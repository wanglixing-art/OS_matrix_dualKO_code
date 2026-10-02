cat("--- rtools compile test ---\n")
td <- file.path(tempdir(), "cpptest")
dir.create(td, showWarnings = FALSE)
cpp <- file.path(td, "hello.cpp")
writeLines('extern "C" { int hello_cpp(){ return 42; } }', cpp)
r <- system2(file.path(R.home("bin"), "Rcmd.exe"),
             args = c("SHLIB", cpp),
             stdout = TRUE, stderr = TRUE)
cat(paste(r, collapse = "\n"), "\n")
dll <- file.path(td, paste0("hello", .Platform$dynlib.ext))
cat("dll exists:", file.exists(dll), "\n")
if (file.exists(dll)) {
  dyn.load(dll)
  cat("call result:", .C("hello_cpp")[[1]], "\n")
}
cat("--- done ---\n")
