# Save the C++ code to a file
cpp_file <- tempfile(tmpdir = ".", fileext = ".cpp")
writeLines(readLines("example.cpp"), cpp_file)

deps_paths <- function(cxxargs, package) {
  # code from inline package
  dir <- system.file("include", package = package)
  if (.Platform$OS.type == "windows") {
    dir <- utils::shortPathName(normalizePath(dir))
  }
  cxxargs <- c(paste("-I", dir, sep = ""), cxxargs)
}
cxxargs <- "-std=c++20"
cxxargs <- deps_paths(cxxargs, "Rcpp")
cxxargs <- deps_paths(cxxargs, "RcppArmadillo")
cxxargs <- deps_paths(cxxargs, "ast2ast")
args <- paste(cxxargs, collapse = " ")

Sys.setenv(PKG_CXXFLAGS = args)


# Compile the shared library
output <- system2("R", args = c("CMD", "SHLIB", cpp_file), stdout = TRUE, stderr = TRUE)
cat(output, sep = "\n")


# Determine the name of the shared library
lib_name <- function(f) {
  if (.Platform$OS.type == "windows") {
    return(paste0(f, ".dll"))
  }
  return(paste0(f, ".so"))
}

# Remove file extension and construct library path
f_no_ext <- sub("\\.[^.]*$", "", basename(cpp_file))
lib <- lib_name(f_no_ext)
lib <- normalizePath(lib)

# Load the shared library
dyn.load(lib)

# Call the function
.Call("test")

# Clean up
fs <- list.files()
unlink(fs[grepl("\\.cpp", fs)])
unlink(fs[grepl("\\.o", fs)])
unlink(fs[grepl("\\.so", fs)])
