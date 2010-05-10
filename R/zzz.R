.onLoad <- function(lib, pkg) {
	pkgVersion <- packageDescription(pkg)$Version
	msg <- paste("\nWelcome to", pkg, "version", pkgVersion, "\n")
	packageStartupMessage(msg)
	.initial()
}
