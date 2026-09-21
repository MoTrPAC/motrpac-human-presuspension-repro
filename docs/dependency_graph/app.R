# motrpac-human-presuspension-repro dependency-graph explorer (Shiny)
# -----------------------------------------------------------------------------
# Interactive companion to build_depgraph.R. Reads the node/edge tables that
# build_depgraph.R emits (depgraph_nodes.csv / depgraph_edges.csv) and lets you
# pick any node and trace it in either direction:
#   UPSTREAM   -- everything that has to be correct for the node to reproduce
#   DOWNSTREAM -- everything that consumes the node (e.g. which freeze files a
#                 shared helper reaches, or which bucket a stem's output lands in)
#   BOTH       -- upstream on the left, the node in the middle, downstream on the
#                 right, in one dataflow picture.
#
# Run:
#   shiny::runApp("docs/dependency_graph")
# or from inside this folder:
#   shiny::runApp()
#
# Needs: shiny, visNetwork, igraph, dplyr, DT. It reads the two CSVs only — no
# pipeline run, no gitignored build outputs, no consortium data.
# -----------------------------------------------------------------------------

library(shiny)
library(visNetwork)
library(igraph)
library(dplyr)
library(DT)

# --- locate the data next to this script ------------------------------------
app_dir = tryCatch(dirname(knitr::current_input(dir = TRUE)), error = function(e) NULL)
if (is.null(app_dir) || !nzchar(app_dir)) {
  app_dir = getwd()
}
nodes_path = file.path(app_dir, "depgraph_nodes.csv")
edges_path = file.path(app_dir, "depgraph_edges.csv")
if (!file.exists(nodes_path)) nodes_path = "depgraph_nodes.csv"
if (!file.exists(edges_path)) edges_path = "depgraph_edges.csv"
if (!file.exists(nodes_path))
  stop("depgraph_nodes.csv not found — run: Rscript docs/dependency_graph/build_depgraph.R")

nodes_raw = read.csv(nodes_path, stringsAsFactors = FALSE)
edges_raw = read.csv(edges_path, stringsAsFactors = FALSE)

# keep only edges whose endpoints are real nodes
valid_ids = nodes_raw$id
edges_raw = edges_raw %>%
  dplyr::filter(from %in% valid_ids, to %in% valid_ids)

# --- canonical DEPENDENCY orientation ---------------------------------------
# The CSVs are written in DATAFLOW orientation (the arrow points at whatever
# comes next). A dependency edge a -> b means "a depends on b":
#   produces    (stem -> freeze)   : the file depends on the stem      -> reverse
#   consumed_by (freeze -> step)   : the step depends on the file      -> reverse
#   feeds       (obj -> obj)       : the later object depends on it    -> reverse
#   runs        (driver -> step)   : the step depends on its driver    -> reverse
#   orders      (target -> target) : the later target depends on it    -> reverse
#   requires    (input -> gate)    : the gate depends on the input     -> reverse
#   defines     (file -> fn)       : the fn lives in the file          -> reverse
#   sourced_by  (lib -> script)    : the script depends on the lib     -> reverse
#   checked_by  (asset -> preflight): preflight depends on the asset   -> reverse
#   calls       (caller -> callee) : the caller depends on the callee  -> keep
# UPSTREAM of X = follow dependency edges from X (mode = "out").
# DOWNSTREAM of X = who depends on X                (mode = "in").
dep_edges = edges_raw %>%
  dplyr::mutate(
    dep_from = ifelse(type == "calls", from, to),
    dep_to   = ifelse(type == "calls", to,   from)
  ) %>%
  dplyr::select(dep_from, dep_to, type)

g_dep = igraph::graph_from_data_frame(
  d = dep_edges %>% dplyr::select(dep_from, dep_to),
  vertices = data.frame(name = nodes_raw$id, stringsAsFactors = FALSE),
  directed = TRUE
)

# --- display styling (mirrors build_depgraph.R / README legend) -------------
stage_color = c(
  orchestration = "#455A64",   # slate  — Makefile
  stage0        = "#00897B",   # teal   — preflight
  stage1        = "#7E57C2",   # purple — build data
  stage2        = "#FB8C00",   # orange — upload
  stage3        = "#C2185B",   # pink   — downstream packages
  external      = "#9E9E9E"    # grey   — outside this repo
)
type_shape = c(
  make_target = "database", stage_driver = "box", step = "square",
  stem = "square", test_script = "star", lib_file = "square",
  lib_fn = "triangleDown", data_object = "dot", freeze_file = "dot",
  preflight_builder = "square", preflight_object = "dot",
  source = "diamond", gate = "hexagon", gate_input = "diamond", sink = "box"
)
type_labels = c(
  make_target = "make target", stage_driver = "stage driver",
  step = "build step", stem = "R generator stem", test_script = "test script",
  lib_file = "shared lib file", lib_fn = "shared lib function",
  data_object = "data object group", freeze_file = "freeze output",
  preflight_builder = "preflight builder", preflight_object = "preflight asset",
  source = "external source", gate = "required-inputs gate",
  gate_input = "required input", sink = "downstream sink"
)
node_shape = function(t) unname(ifelse(is.na(type_shape[t]), "dot", type_shape[t]))
node_color = function(s) unname(ifelse(is.na(stage_color[s]), "#9E9E9E", stage_color[s]))

ome_choices = sort(unique(nodes_raw$ome))
root_type_choices = c(
  "freeze outputs"  = "freeze_file",
  "R generator stems" = "stem",
  "build steps"     = "step",
  "data objects"    = "data_object",
  "shared lib fns"  = "lib_fn",
  "external sources" = "source",
  "make targets"    = "make_target",
  "all node types"  = "all"
)

# --- UI ----------------------------------------------------------------------
ui = fluidPage(
  titlePanel("motrpac-human-presuspension-repro dependency explorer"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      radioButtons("direction", "Direction",
                   choices = c("Both (upstream + downstream)" = "both",
                               "Upstream (what it depends on)" = "up",
                               "Downstream (what depends on it)" = "down"),
                   selected = "up"),
      selectInput("root_type", "Root list: node type",
                  choices = root_type_choices, selected = "freeze_file"),
      selectInput("ome", "Filter root list by ome",
                  choices = c("all", ome_choices), selected = "all"),
      selectizeInput("root", "Focal node", choices = NULL,
                     options = list(placeholder = "select a node...")),
      sliderInput("hops", "Max hops", min = 1, max = 14, value = 5, step = 1),
      checkboxInput("show_calls", "Include shared-lib calls / defines",
                    value = FALSE),
      checkboxInput("show_control", "Include make/driver control flow",
                    value = TRUE),
      hr(),
      helpText("Click a node in the graph to re-root on it.",
               "Color = stage, shape = type. Left = upstream, right = downstream."),
      htmlOutput("stats")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel("Graph", visNetworkOutput("net", height = "720px")),
        tabPanel("Dependency table", DT::dataTableOutput("tbl"))
      )
    )
  )
)

# --- server ------------------------------------------------------------------
server = function(input, output, session) {

  # populate the focal-node dropdown based on node-type + ome filter
  root_choices = reactive({
    df = nodes_raw
    if (input$root_type != "all") df = df %>% dplyr::filter(type == input$root_type)
    if (input$ome != "all")       df = df %>% dplyr::filter(ome == input$ome)
    df = df %>% dplyr::arrange(label)
    stats::setNames(df$id, sprintf("%s  [%s]", df$label, df$stage))
  })

  observeEvent(root_choices(), {
    ch = root_choices()
    sel = isolate(input$root)
    if (is.null(sel) || !(sel %in% ch)) sel = if (length(ch)) ch[[1]] else character(0)
    updateSelectizeInput(session, "root", choices = ch, selected = sel, server = TRUE)
  }, ignoreNULL = FALSE)

  # reachable set (with signed hop distance + relation) from the focal node
  reach = reactive({
    req(input$root)
    root = input$root
    if (!(root %in% igraph::V(g_dep)$name)) return(NULL)

    grab = function(mode) {
      d = igraph::distances(g_dep, v = root, mode = mode)[1, ]
      d = d[is.finite(d) & d <= input$hops]
      d = d[names(d) != root]
      if (!length(d)) return(NULL)
      data.frame(id = names(d), dist = as.integer(d), stringsAsFactors = FALSE)
    }
    up   = if (input$direction %in% c("up", "both"))   grab("out") else NULL
    down = if (input$direction %in% c("down", "both")) grab("in")  else NULL
    if (!is.null(up))   up$rel   = "upstream"
    if (!is.null(down)) down$rel = "downstream"

    focal = data.frame(id = root, dist = 0L, rel = "focal", stringsAsFactors = FALSE)
    out = dplyr::bind_rows(focal, up, down)
    # a node reachable both ways (rare, only via cycles): keep the closer hop
    out = out %>%
      dplyr::group_by(id) %>%
      dplyr::slice_min(order_by = dist, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup()
    out
  })

  # build the visNetwork data frames for the reachable subgraph
  vis_data = reactive({
    rr = reach()
    if (is.null(rr) || nrow(rr) == 0) return(NULL)

    nd = nodes_raw %>% dplyr::filter(id %in% rr$id)
    rr = rr %>% dplyr::filter(id %in% nd$id)

    # layout levels: upstream to the left, focal middle, downstream right
    maxup = max(rr$dist[rr$rel == "upstream"], 0)
    lev_of = function(rel, dist) {
      ifelse(rel == "upstream", maxup - dist,
             ifelse(rel == "focal", maxup, maxup + dist))
    }

    vn = nd %>%
      dplyr::left_join(rr, by = "id") %>%
      dplyr::mutate(
        shape = node_shape(type),
        color = node_color(stage),
        level = lev_of(rel, dist),
        borderWidth = ifelse(id == input$root, 4, ifelse(hotspot, 3, 1)),
        size = ifelse(id == input$root, 34, 14 + pmin(n_downstream, 30) * 0.6),
        font.size = ifelse(id == input$root, 26, 16),
        title = sprintf(
          "<b>%s</b><br>%s | %s | %s<br>%s, %d hop(s) | downstream reach: %d%s%s<br><code>%s</code>",
          label, type_labels[type], stage, ome, rel, dist, n_downstream,
          ifelse(hotspot, " | <b>hotspot</b>", ""),
          ifelse(nzchar(detail), paste0("<br>", detail), ""), id
        )
      ) %>%
      dplyr::rename(group = stage) %>%
      dplyr::select(id, label, shape, color, level, borderWidth, size,
                    font.size, title, group)
    vn$color[vn$id == input$root] = "#E53935"  # focal = red halo

    keep = nd$id
    # subgraph edges in dataflow orientation (arrows point toward the consumer),
    # so the picture reads sources -> focal -> bucket left to right
    ke = edges_raw %>% dplyr::filter(from %in% keep, to %in% keep)
    edge_col = c(orders = "#455A64", runs = "#5C6BC0", produces = "#2E7D32",
                 consumed_by = "#546E7A", feeds = "#EF6C00", requires = "#C62828",
                 tested_by = "#00838F", checked_by = "#8D6E63",
                 calls = "#B0BEC5", defines = "#CFD8DC", sourced_by = "#90A4AE")
    if (!isTRUE(input$show_calls))
      ke = ke %>% dplyr::filter(!type %in% c("calls", "defines"))
    if (!isTRUE(input$show_control))
      ke = ke %>% dplyr::filter(!type %in% c("runs", "orders"))
    ve = ke %>%
      dplyr::transmute(
        from = from, to = to, arrows = "to",
        color = unname(ifelse(is.na(edge_col[type]), "#B0BEC5", edge_col[type])),
        dashes = type %in% c("calls", "defines", "sourced_by"),
        title = type
      )
    list(nodes = vn, edges = ve)
  })

  output$net = renderVisNetwork({
    vd = vis_data()
    validate(need(!is.null(vd) && nrow(vd$nodes) > 1,
                  "Select a focal node that has dependencies in this direction."))
    visNetwork(vd$nodes, vd$edges) %>%
      visHierarchicalLayout(direction = "LR", sortMethod = "directed",
                            levelSeparation = 220, nodeSpacing = 120) %>%
      visNodes(borderWidthSelected = 5) %>%
      visEdges(smooth = list(enabled = TRUE, type = "cubicBezier")) %>%
      visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
                 nodesIdSelection = FALSE) %>%
      visInteraction(hover = TRUE, tooltipDelay = 120) %>%
      visEvents(click = "function(p){
        if (p.nodes.length > 0) { Shiny.setInputValue('clicked_node', p.nodes[0], {priority:'event'}); }
      }")
  })

  # click a node -> re-root on it (open up the type filter so any node is selectable)
  observeEvent(input$clicked_node, {
    nid = input$clicked_node
    if (!(nid %in% nodes_raw$id)) return()
    ntype = nodes_raw$type[nodes_raw$id == nid]
    if (input$root_type != "all" && input$root_type != ntype) {
      updateSelectInput(session, "root_type", selected = "all")
    }
    if (input$ome != "all" && nodes_raw$ome[nodes_raw$id == nid] != input$ome) {
      updateSelectInput(session, "ome", selected = "all")
    }
    updateSelectizeInput(session, "root", selected = nid, server = TRUE)
  })

  output$tbl = DT::renderDataTable({
    rr = reach(); req(rr)
    rr %>%
      dplyr::filter(rel != "focal") %>%
      dplyr::inner_join(nodes_raw, by = "id") %>%
      dplyr::transmute(
        relation = rel,
        hops = ifelse(rel == "upstream", -dist, dist),
        node = label, type = type_labels[type], stage = stage, ome = ome,
        downstream_reach = n_downstream, hotspot = hotspot, id = id
      ) %>%
      dplyr::arrange(relation, abs(hops), dplyr::desc(downstream_reach))
  }, options = list(pageLength = 25, order = list()), rownames = FALSE)

  output$stats = renderUI({
    rr = reach(); req(rr)
    root = input$root
    focal_lab = nodes_raw$label[nodes_raw$id == root]
    n_up   = sum(rr$rel == "upstream")
    n_down = sum(rr$rel == "downstream")
    frz_down = nodes_raw %>%
      dplyr::filter(id %in% rr$id[rr$rel == "downstream"], type == "freeze_file") %>%
      nrow()
    gated_up = nodes_raw %>%
      dplyr::filter(id %in% rr$id[rr$rel == "upstream"], origin == "gated") %>%
      nrow()
    HTML(paste0(
      "<b>", focal_lab, "</b><br>within ", input$hops, " hops:<br>",
      "&nbsp;&nbsp;upstream deps: <b>", n_up, "</b><br>",
      "&nbsp;&nbsp;downstream consumers: <b>", n_down, "</b><br>",
      "&nbsp;&nbsp;&#9642; freeze outputs downstream: <b>", frz_down, "</b><br>",
      "&nbsp;&nbsp;&#9642; consortium-gated inputs upstream: <b>", gated_up, "</b>"
    ))
  })
}

shinyApp(ui, server)
