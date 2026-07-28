\name{bp_check_versions}
\alias{bp_check_versions}
\title{Check Audited Package Versions}
\usage{
bp_check_versions(error = FALSE)
}
\arguments{
\item{error}{Logical scalar. If \code{TRUE}, throw an error when the installed
versions do not match a complete audited stack. If \code{FALSE}, return the
check table without error.}
}
\value{
A data frame with package names, installed versions, expected versions,
individual version support, complete-stack match status, and the matched stack
name. A version can be individually supported while \code{ok} is \code{FALSE}
when the installed versions form an unaudited mixture.
}
\description{
Reports whether the installed R, edgeR, limma, and BPCells versions match one
complete stack audited for exact streamed edgeR compatibility.
}
