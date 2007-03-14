\name{samplesMDS}
\alias{samplesMDS}
\alias{speciesMDS}
\title{Create an MDS plot for a RiboSort object}
\description{A object produced by the function RiboSort is a matrix with columns representing samples and rows representing species. \code{samplesMDS} and \code{speciesMDS} produce MDS plots of samples and species respectively.}

\usage{
   samplesMDS(x, dissimilarity = "euclidean", type = "non-metric", labels=TRUE)
   speciesMDS(x, dissimilarity = "euclidean", type = "non-metric", labels=TRUE)
}
\arguments{
  \item{x}{A numeric matrix, dataframe or a "RiboSort" object, ie. an object produced by the RiboSort function}
  \item{dissimilarity}{The distance measure to be used computing dissimilarities between samples or species. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". Any unambiguous substring can be given.}
  \item{type}{The type of Multi-dimensional Scaling to be used. }
  \item{labels}{A logical argument. When TRUE, labels are included in the MDS plot. If FALSE, labels are omitted.}
}
\details{
  These functions allow for simple and quick production of Multi-dimensional Scaling plots. A choice of dissimilarities is available, aswell as the option to specify the desired type of multi-dimensional scaling. Note that at least three samples are required to produce a two-dimensional plot.
  
When deciding upon a dissimilarity measure, consult the description of available distance measures in the helpfile \code{?dist}.

When choosing which type of multi-dimensional scaling to use, the following elaboration on the options available may aid your choice.
\item \code{classical :}  Classical multi-dimensional scaling of a data matrix. Also known as principal coordinates analysis.
\item \code{sammon :}  Sammon's non-linear mapping, also known as metric least squares  multi-dimensional scaling.
\item \code{non-metric :}  Kruskal's form of non-metric multidimensional scaling.
}
\value{
  If \code{type = "classical"}, a matrix with two columns whose rows give the coordinates of the points chosen to represent the dissimilarities.


  If \code{type = "sammon"} or \code{"non-metric"}, a list containing the following two components:
  \item{points}{A two-column vector of the fitted configuration.}
  \item{stress}{The final stress achieved (in percent).}
   
}
\author{Una Scallan, \email{una.scallan@ucd.ie}}
\note{
   For a more interactive and less limited approach to producing multi-dimensional scaling graphs, see the functions listed below under: See Also.
}
\section{Warning}{
   If there appears an error stating that the initial 'configuration has duplicates', this implies that the dissimilarity specified is not suitable. Its use in producing a distance matrix, led to the creation of duplicates in the initial configuration (computed via Classical MDS), and thus prevented an MDS plot from being created. Try alternative dissimilarity measures.

   If there appears an error stating that the k in 'cmdscale(d, k) must be a value in {1,2,...,n-1}', then this implies that you have less than three samples represented in the object (usually the RiboSort object, x) supplied to samplesMDS. It is not sensible to create a two-dimensional plot when less than three samples are being compared. The difference between two samples can easily be illustrated on a one-dimensional graph (i.e. a line).
}

\seealso{ \code{\link{RiboSort}}, \code{link{dist}}, \code{\link{cmdscale}}, \code{\link[MASS]{sammon}},\ code{\link[MASS]{isoMDS}} }

\examples{

#---Example 1
data(demo.beckman)
file.copy(system.file(demo.beckman, package = "RiboSort"), ".") #copying demo datafiles provided with RiboSort into the current working directory
y=RiboSort(data=demo.beckman, dataformat="beckman", dye="D4")
samplesMDS(y,dissimilarity="euclidean", type="non-metric", labels=TRUE)

#---Example 2
data(demo.abi.multiple)
file.copy(system.file(demo.abi.multiple, package = "RiboSort"), ".")
y=RiboSort(data=demo.abi.multiple, dataformat="abimultiple")
samplesMDS(y,dissimilarity="manhattan", type="sammon", labels=TRUE)

#---Example 3
data(demo.standard)
file.copy(system.file(demo.standard, package = "RiboSort"), ".")
y=RiboSort(data=demo.standard, dataformat="standard", repeats=2, mergerepeats="presentinone")
par(mfrow=c(1,2))	 #setting graphics window to display 2 graphs
samplesMDS(y,dissimilarity="canberra", type="classical", labels=TRUE)
speciesMDS(y,dissimilarity="canberra", type="classical", labels=TRUE)		

}
\keyword{multivariate}
