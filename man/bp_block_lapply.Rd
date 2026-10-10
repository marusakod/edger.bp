\name{bp_block_lapply}
\alias{bp_block_lapply}
\title{Apply a block callback using the active worker pool}
\usage{bp_block_lapply(x, workers, FUN, context = NULL)}
\arguments{
\item{x}{List or vector of tasks.}
\item{workers}{Number of workers, matching an active pool's size.}
\item{FUN}{Function of one task. Captured disk-backed matrices must equal
\code{pool$counts}. Do not launch nested parallel jobs inside the callback.}
\item{context}{Optional label for error and progress messages.}
}
\description{Uses the current PSOCK scope, or the usual numeric-worker backend
outside a scope. Only required callback bindings are sent to socket workers.}
\value{A list of results in input order. Task errors include the task and worker identity.}
