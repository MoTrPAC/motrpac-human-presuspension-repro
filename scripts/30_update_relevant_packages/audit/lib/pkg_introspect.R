# Static introspection helpers for an R source package.
#
# Everything here reads source text and DESCRIPTION/NAMESPACE metadata. Nothing
# installs, loads or attaches the package under audit, so the audit runs on a
# package that does not currently install — which is the case an audit is for.
#
# Function definitions are the one exception: `eval`ing a `name <- function(...)`
# expression binds the closure without running its body, and that is what
# codetools::findGlobals needs to resolve call edges. No other top-level code is
# evaluated.

# ---- Package layout --------------------------------------------------------

pkg_r_files <- function(pkg_root) {
  r_dir <- file.path(pkg_root, "R")
  if (!dir.exists(r_dir)) return(character(0))
  sort(list.files(r_dir, pattern = "\\.[RrSsQq]$", full.names = TRUE))
}

pkg_data_objects <- function(pkg_root) {
  data_dir <- file.path(pkg_root, "data")
  if (!dir.exists(data_dir)) return(character(0))
  files <- list.files(data_dir, pattern = "\\.(rda|RData|rds)$", full.names = FALSE)
  sort(unique(sub("\\.(rda|RData|rds)$", "", files)))
}

pkg_man_topics <- function(pkg_root) {
  man_dir <- file.path(pkg_root, "man")
  if (!dir.exists(man_dir)) return(character(0))
  sort(list.files(man_dir, pattern = "\\.Rd$", full.names = FALSE))
}

# The \alias{} entries of every man/*.Rd, which is what `?name` actually resolves
# against. A data object is documented iff its name appears as an alias.
pkg_man_aliases <- function(pkg_root) {
  man_dir <- file.path(pkg_root, "man")
  if (!dir.exists(man_dir)) return(data.frame(rd = character(0), alias = character(0)))
  rds <- list.files(man_dir, pattern = "\\.Rd$", full.names = TRUE)
  out <- lapply(rds, function(rd) {
    txt <- readLines(rd, warn = FALSE)
    hits <- regmatches(txt, gregexpr("\\\\alias\\{[^}]*\\}", txt))
    alias <- gsub("^\\\\alias\\{|\\}$", "", unlist(hits))
    if (!length(alias)) return(NULL)
    data.frame(rd = basename(rd), alias = alias, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out)
  if (is.null(out)) data.frame(rd = character(0), alias = character(0)) else out
}

# Rd topics carrying \docType{data} — the documentation form R CMD check requires
# for anything shipped in data/.
#
# Returns aliases and topic names separately. Under @rdname, one Rd documents many
# objects and carries an alias for each PLUS one for the topic itself; that topic
# alias (e.g. QC_RESULTS covering 42 matrices) names no object in data/ and must
# not be mistaken for documentation of a deleted one.
pkg_data_doc_aliases <- function(pkg_root) {
  man_dir <- file.path(pkg_root, "man")
  empty <- list(aliases = character(0), topics = character(0))
  if (!dir.exists(man_dir)) return(empty)
  rds <- list.files(man_dir, pattern = "\\.Rd$", full.names = TRUE)
  aliases <- character(0)
  topics <- character(0)
  for (rd in rds) {
    txt <- paste(readLines(rd, warn = FALSE), collapse = "\n")
    if (!grepl("\\\\docType\\{data\\}", txt)) next
    hits <- regmatches(txt, gregexpr("\\\\alias\\{[^}]*\\}", txt))
    aliases <- c(aliases, gsub("^\\\\alias\\{|\\}$", "", unlist(hits)))
    nm <- regmatches(txt, regexpr("\\\\name\\{[^}]*\\}", txt))
    if (length(nm)) topics <- c(topics, gsub("^\\\\name\\{|\\}$", "", nm))
  }
  list(aliases = sort(unique(aliases)), topics = sort(unique(topics)))
}

# ---- DESCRIPTION -----------------------------------------------------------

read_description <- function(pkg_root) {
  path <- file.path(pkg_root, "DESCRIPTION")
  if (!file.exists(path)) stop("no DESCRIPTION at ", pkg_root)
  as.data.frame(read.dcf(path), stringsAsFactors = FALSE)
}

# Split one dependency field into bare package names, dropping version
# constraints and the implicit R entry.
description_deps <- function(desc, field = c("Depends", "Imports", "Suggests",
                                             "LinkingTo", "Enhances")) {
  field <- match.arg(field)
  if (!field %in% names(desc)) return(character(0))
  raw <- desc[[field]]
  if (is.na(raw) || !nzchar(raw)) return(character(0))
  parts <- trimws(strsplit(raw, ",")[[1]])
  parts <- sub("\\s*\\(.*\\)$", "", parts)
  parts <- trimws(parts)
  parts <- parts[nzchar(parts) & parts != "R"]
  sort(unique(parts))
}

# ---- NAMESPACE -------------------------------------------------------------

parse_namespace <- function(pkg_root) {
  path <- file.path(pkg_root, "NAMESPACE")
  empty <- list(exports = character(0), imports = character(0),
                import_from = data.frame(pkg = character(0), fn = character(0)),
                generated = FALSE)
  if (!file.exists(path)) return(empty)
  txt <- readLines(path, warn = FALSE)

  grab <- function(pattern) {
    hits <- regmatches(txt, regexpr(pattern, txt))
    hits[nzchar(hits)]
  }
  exports <- gsub('^export\\(|\\)$|"', "", grab('^export\\([^)]*\\)'))
  s3 <- gsub('^S3method\\(|\\)$|"', "", grab('^S3method\\([^)]*\\)'))
  imports <- gsub('^import\\(|\\)$|"', "", grab('^import\\([^)]*\\)'))

  from_lines <- grab('^importFrom\\([^)]*\\)')
  from <- gsub('^importFrom\\(|\\)$|"', "", from_lines)
  from_split <- strsplit(from, "\\s*,\\s*")
  import_from <- do.call(rbind, lapply(from_split, function(x) {
    if (length(x) < 2) return(NULL)
    data.frame(pkg = x[1], fn = x[-1], stringsAsFactors = FALSE)
  }))
  if (is.null(import_from)) {
    import_from <- data.frame(pkg = character(0), fn = character(0))
  }

  list(exports = sort(unique(trimws(exports))),
       s3methods = sort(unique(trimws(s3))),
       imports = sort(unique(trimws(imports))),
       import_from = import_from,
       generated = any(grepl("Generated by roxygen2", txt, fixed = TRUE)))
}

# ---- Function definitions and the call graph -------------------------------

# Top-level `name <- function(...)` definitions, one row per definition, with the
# file and line it came from. Duplicated names across files are kept: a name
# defined twice is itself a finding.
pkg_function_defs <- function(pkg_root) {
  files <- pkg_r_files(pkg_root)
  rows <- lapply(files, function(f) {
    exprs <- tryCatch(parse(f, keep.source = TRUE),
                      error = function(e) NULL)
    if (is.null(exprs)) {
      warning("could not parse ", f, call. = FALSE)
      return(NULL)
    }
    refs <- utils::getSrcref(exprs)
    out <- lapply(seq_along(exprs), function(i) {
      e <- exprs[[i]]
      if (!is.call(e) || length(e) < 3L) return(NULL)
      op <- as.character(e[[1]])
      if (!op %in% c("<-", "=", "<<-")) return(NULL)
      if (!is.name(e[[2]])) return(NULL)
      rhs <- e[[3]]
      if (!(is.call(rhs) && identical(as.character(rhs[[1]]), "function"))) return(NULL)
      data.frame(name = as.character(e[[2]]),
                 file = f,
                 line = if (!is.null(refs) && length(refs) >= i) refs[[i]][1L] else NA_integer_,
                 stringsAsFactors = FALSE)
    })
    do.call(rbind, out)
  })
  rows <- do.call(rbind, rows)
  if (is.null(rows)) {
    return(data.frame(name = character(0), file = character(0), line = integer(0)))
  }
  rows$file <- normalizePath(rows$file, mustWork = FALSE)
  rows[order(rows$name), ]
}

# Bind every top-level function literal into one environment so findGlobals can
# resolve its free variables. The body is never run.
pkg_function_env <- function(pkg_root) {
  env <- new.env(parent = baseenv())
  for (f in pkg_r_files(pkg_root)) {
    exprs <- tryCatch(parse(f, keep.source = FALSE), error = function(e) NULL)
    if (is.null(exprs)) next
    for (e in exprs) {
      if (!is.call(e) || length(e) < 3L) next
      if (!as.character(e[[1]]) %in% c("<-", "=", "<<-")) next
      if (!is.name(e[[2]])) next
      rhs <- e[[3]]
      if (!(is.call(rhs) && identical(as.character(rhs[[1]]), "function"))) next
      tryCatch(assign(as.character(e[[2]]), eval(rhs, envir = env), envir = env),
               error = function(err) NULL)
    }
  }
  env
}

# Edges caller -> callee, restricted to functions this package defines. Calls out
# to other packages are dropped here; declared_but_unused() handles those.
#
# Both halves of findGlobals are used. A helper passed by name rather than called
# — `heatmap_color_fun = .feature_color_function`, or the value argument of
# assignInNamespace — lands in $variables, not $functions. Reading only
# $functions reports every callback and every injected function as unreferenced.
pkg_call_graph <- function(pkg_root) {
  env <- pkg_function_env(pkg_root)
  defined <- ls(env, all.names = TRUE)
  if (!length(defined)) {
    return(data.frame(caller = character(0), callee = character(0)))
  }
  rows <- lapply(defined, function(nm) {
    fn <- get(nm, envir = env)
    if (!is.function(fn)) return(NULL)
    g <- tryCatch(codetools::findGlobals(fn, merge = FALSE),
                  error = function(e) list(functions = character(0),
                                           variables = character(0)))
    called <- unique(c(g$functions, g$variables))
    called <- intersect(called, defined)
    called <- setdiff(called, nm)          # self-recursion tracked separately
    if (!length(called)) return(NULL)
    data.frame(caller = nm, callee = called, stringsAsFactors = FALSE)
  })
  rows <- do.call(rbind, rows)
  if (is.null(rows)) data.frame(caller = character(0), callee = character(0)) else rows
}

pkg_self_recursive <- function(pkg_root) {
  env <- pkg_function_env(pkg_root)
  defined <- ls(env, all.names = TRUE)
  hits <- vapply(defined, function(nm) {
    fn <- get(nm, envir = env)
    if (!is.function(fn)) return(FALSE)
    called <- tryCatch(codetools::findGlobals(fn, merge = FALSE)$functions,
                       error = function(e) character(0))
    nm %in% called
  }, logical(1))
  sort(defined[hits])
}

# Every strongly connected component with more than one node, i.e. every genuine
# call cycle. Tarjan, iterative, so a deep graph cannot blow the R stack.
find_cycles <- function(edges) {
  nodes <- sort(unique(c(edges$caller, edges$callee)))
  if (!length(nodes)) return(list())
  adj <- split(edges$callee, factor(edges$caller, levels = nodes))
  adj <- lapply(adj, function(x) if (is.null(x)) character(0) else unique(x))

  index <- setNames(rep(NA_integer_, length(nodes)), nodes)
  lowlink <- setNames(rep(NA_integer_, length(nodes)), nodes)
  on_stack <- setNames(rep(FALSE, length(nodes)), nodes)
  stack <- character(0)
  counter <- 0L
  components <- list()

  for (root in nodes) {
    if (!is.na(index[[root]])) next
    # work holds (node, next-child-position) frames
    work <- list(list(v = root, i = 1L))
    while (length(work)) {
      frame <- work[[length(work)]]
      v <- frame$v
      if (frame$i == 1L) {
        counter <- counter + 1L
        index[[v]] <- counter
        lowlink[[v]] <- counter
        stack <- c(stack, v)
        on_stack[[v]] <- TRUE
      }
      children <- adj[[v]]
      recursed <- FALSE
      while (frame$i <= length(children)) {
        w <- children[[frame$i]]
        frame$i <- frame$i + 1L
        work[[length(work)]] <- frame
        if (is.na(index[[w]])) {
          work[[length(work) + 1L]] <- list(v = w, i = 1L)
          recursed <- TRUE
          break
        } else if (on_stack[[w]]) {
          lowlink[[v]] <- min(lowlink[[v]], index[[w]])
        }
      }
      if (recursed) next
      work[[length(work)]] <- frame
      if (frame$i > length(children)) {
        if (lowlink[[v]] == index[[v]]) {
          pos <- which(stack == v)
          pos <- pos[length(pos)]
          comp <- stack[pos:length(stack)]
          stack <- if (pos > 1L) stack[seq_len(pos - 1L)] else character(0)
          on_stack[comp] <- FALSE
          if (length(comp) > 1L) components[[length(components) + 1L]] <- sort(comp)
        }
        work[[length(work)]] <- NULL
        if (length(work)) {
          parent <- work[[length(work)]]
          lowlink[[parent$v]] <- min(lowlink[[parent$v]], lowlink[[v]])
        }
      }
    }
  }
  components
}

# ---- Namespace-qualified calls ---------------------------------------------

# Every `pkg` appearing as `pkg::fn` or `pkg:::fn` in R/, taken from the parser's
# token stream rather than by grepping the text. R tags that operand
# SYMBOL_PACKAGE, so comments, roxygen and string literals cannot produce a hit —
# grepping for "pkg::" finds all three and invents dependencies that do not exist.
pkg_qualified_calls <- function(pkg_root) {
  files <- pkg_r_files(pkg_root)
  rows <- lapply(files, function(f) {
    pd <- tryCatch({
      exprs <- parse(f, keep.source = TRUE)
      utils::getParseData(exprs)
    }, error = function(e) NULL)
    if (is.null(pd) || !nrow(pd)) return(NULL)
    sel <- pd$token == "SYMBOL_PACKAGE"
    if (!any(sel)) return(NULL)
    data.frame(pkg = gsub('^"|"$|^`|`$', "", pd$text[sel]),
               file = f, line = pd$line1[sel], stringsAsFactors = FALSE)
  })
  rows <- do.call(rbind, rows)
  if (is.null(rows)) {
    data.frame(pkg = character(0), file = character(0), line = integer(0))
  } else {
    rows
  }
}

# Every name used as a SYMBOL in R/ — i.e. actually referenced by code.
#
# String constants are excluded by construction, which is the point: a roxygen
# data block ends with the bare string "OBJECT_NAME" naming the object it
# documents. Grepping R/ for a data object's name therefore always finds at least
# that stub and reports every shipped object as used. Documentation is not usage.
pkg_symbol_uses <- function(pkg_root) {
  names(pkg_symbol_use_counts(pkg_root))
}

# How many times each name is referenced in R/, counting SYMBOL tokens and
# excluding the line each function is defined on.
#
# This is the reliable usage signal, and it is why the call graph is not used for
# it. codetools classifies a name as local the moment it is assigned, so a helper
# that is read and then rebound — `environment(.custom_block_set_env) <- ns` —
# disappears from findGlobals entirely and reads as unreferenced. Counting tokens
# and subtracting the definition has no such blind spot.
pkg_symbol_use_counts <- function(pkg_root, def_lines = NULL) {
  if (is.null(def_lines)) {
    defs <- pkg_function_defs(pkg_root)
    def_lines <- paste(defs$file, defs$line, sep = ":")
  }
  files <- pkg_r_files(pkg_root)
  out <- lapply(files, function(f) {
    pd <- tryCatch({
      exprs <- parse(f, keep.source = TRUE)
      utils::getParseData(exprs)
    }, error = function(e) NULL)
    if (is.null(pd) || !nrow(pd)) return(NULL)
    sel <- pd$token %in% c("SYMBOL", "SYMBOL_FUNCTION_CALL")
    used <- character(0)
    if (any(sel)) {
      key <- paste(normalizePath(f, mustWork = FALSE), pd$line1[sel], sep = ":")
      used <- pd$text[sel][!key %in% def_lines]
    }
    # roxygen inline code: `@format `r .describe_data_format()``. The parser sees
    # a comment, but roxygen evaluates it at document() time, so the call is real.
    # Without this every doc-generating helper looks unused.
    txt <- tryCatch(readLines(f, warn = FALSE), error = function(e) character(0))
    roxy <- grep("^\\s*#'", txt, value = TRUE)
    inline <- unlist(regmatches(roxy, gregexpr("`r [^`]*`", roxy)))
    if (length(inline)) {
      used <- c(used, unlist(regmatches(
        inline, gregexpr("[A-Za-z._][A-Za-z0-9._]*", inline))))
    }
    used
  })
  tab <- table(unlist(out))
  if (!length(tab)) return(integer(0))
  setNames(as.integer(tab), names(tab))
}

# Packages named in library()/require()/requireNamespace()/loadNamespace() calls
# in R/, again via the parser: the package name is the first argument, so read the
# token that follows the opening paren.
pkg_attached_calls <- function(pkg_root,
                               attach_fns = c("library", "require",
                                              "requireNamespace", "loadNamespace")) {
  files <- pkg_r_files(pkg_root)
  out <- lapply(files, function(f) {
    pd <- tryCatch({
      exprs <- parse(f, keep.source = TRUE)
      utils::getParseData(exprs)
    }, error = function(e) NULL)
    if (is.null(pd) || !nrow(pd)) return(NULL)
    pd <- pd[pd$terminal, ]
    pd <- pd[order(pd$line1, pd$col1), ]
    idx <- which(pd$token == "SYMBOL_FUNCTION_CALL" & pd$text %in% attach_fns)
    hits <- character(0)
    for (i in idx) {
      # tokens are ( then the package name, as a symbol or a string
      j <- i + 2L
      if (j <= nrow(pd) && pd$token[j] %in% c("SYMBOL", "STR_CONST")) {
        hits <- c(hits, gsub('^"|"$', "", pd$text[j]))
      }
    }
    if (!length(hits)) NULL else unique(hits)
  })
  sort(unique(unlist(out)))
}

# ---- Cross-package textual references --------------------------------------

# Where does `pkg_root` mention `target_pkg`, and in what capacity? The capacity
# is what decides whether a dependency edge is a hard runtime one or only a
# build/test-time one, and therefore whether a cycle actually breaks an install.
pkg_mentions <- function(pkg_root, target_pkg) {
  files <- list.files(pkg_root, recursive = TRUE, full.names = TRUE,
                      all.files = TRUE, no.. = TRUE)
  files <- files[!grepl("(^|/)(\\.git|\\.Rproj\\.user|renv)(/|$)", files)]
  files <- files[!grepl("\\.(rda|RData|rds|png|jpg|jpeg|pdf|gz|zip|xlsx|Rproj)$",
                        files, ignore.case = TRUE)]
  # Skip directories and anything large enough to be a data blob rather than source.
  info <- file.info(files)
  files <- files[!is.na(info$size) & !info$isdir & info$size < 5e6]

  rows <- lapply(files, function(f) {
    txt <- tryCatch(readLines(f, warn = FALSE), error = function(e) NULL)
    if (is.null(txt) || !length(txt)) return(NULL)
    hit <- grep(target_pkg, txt, fixed = TRUE)
    if (!length(hit)) return(NULL)
    rel <- sub(paste0("^", normalizePath(pkg_root, mustWork = FALSE), "/?"), "",
               normalizePath(f, mustWork = FALSE))
    data.frame(file = rel, line = hit, text = trimws(txt[hit]),
               capacity = classify_mention(rel, txt[hit], target_pkg),
               stringsAsFactors = FALSE)
  })
  rows <- do.call(rbind, rows)
  if (is.null(rows)) {
    data.frame(file = character(0), line = integer(0), text = character(0),
               capacity = character(0))
  } else {
    rows
  }
}

# A mention only creates an install-time dependency if it sits in shipped code or
# in DESCRIPTION/NAMESPACE. data-raw/, tests/, vignettes/ and CI are build-time.
# Vectorised over line_text: one file contributes many hit lines at once.
classify_mention <- function(rel_path, line_text, target_pkg) {
  n <- length(line_text)
  if (identical(rel_path, "DESCRIPTION")) return(rep("DESCRIPTION", n))
  if (identical(rel_path, "NAMESPACE")) return(rep("NAMESPACE", n))
  if (grepl("^R/", rel_path)) {
    # A line in R/ is a runtime edge only if it actually calls into the target;
    # a mention in a comment or a roxygen block creates no dependency at all.
    commented <- grepl("^\\s*#", line_text, useBytes = TRUE)
    qualified <-
      grepl(paste0(target_pkg, "::"), line_text, fixed = TRUE, useBytes = TRUE) |
      grepl(paste0("(library|require)\\(\\s*['\"]?", target_pkg),
            line_text, useBytes = TRUE)
    # The package name reached only as a STRING — asNamespace("Pkg"),
    # requireNamespace("Pkg"), get(x, envir = asNamespace("Pkg")). This is a real
    # runtime edge that no dependency resolver and no R CMD check can see, which
    # is exactly why it gets its own capacity rather than being folded into the
    # qualified case: a cycle hiding here looks acyclic in DESCRIPTION.
    via_string <-
      grepl(paste0("(asNamespace|requireNamespace|loadNamespace|getNamespace)\\(\\s*['\"]",
                   target_pkg, "['\"]"), line_text, useBytes = TRUE) |
      grepl(paste0("['\"]", target_pkg, "['\"]"), line_text, useBytes = TRUE)
    ifelse(commented, "R/ comment or string",
           ifelse(qualified, "R/ runtime call",
                  ifelse(via_string, "R/ runtime via string", "R/ comment or string")))
  } else if (grepl("^data-raw/", rel_path)) {
    rep("data-raw (build-time)", n)
  } else if (grepl("^tests/", rel_path)) {
    rep("tests", n)
  } else if (grepl("^vignettes/", rel_path)) {
    rep("vignettes", n)
  } else if (grepl("^man/", rel_path)) {
    rep("man", n)
  } else if (grepl("^\\.github/", rel_path)) {
    rep("CI", n)
  } else if (grepl("^sandbox/", rel_path)) {
    rep("sandbox", n)
  } else {
    rep("other", n)
  }
}

# ---- .Rbuildignore ---------------------------------------------------------

# TRUE if `entry` (a top-level name) is matched by some .Rbuildignore pattern.
# R applies these as case-insensitive Perl regexes against the path relative to
# the package root, unanchored — the same way R CMD build does.
rbuildignore_covers <- function(pkg_root, entry) {
  path <- file.path(pkg_root, ".Rbuildignore")
  if (!file.exists(path)) return(FALSE)
  pats <- readLines(path, warn = FALSE)
  pats <- trimws(pats)
  pats <- pats[nzchar(pats) & !startsWith(pats, "#")]
  any(vapply(pats, function(p) {
    isTRUE(tryCatch(grepl(p, entry, perl = TRUE, ignore.case = TRUE),
                    error = function(e) FALSE))
  }, logical(1)))
}

# ---- Reporting -------------------------------------------------------------

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(path)
}
