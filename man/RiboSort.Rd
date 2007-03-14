\name{RiboSort}
\alias{RiboSort}
\title{Preparation of Microbial Community Profiles for Statistical Analyses}
\description{RiboSort is a function for convenient editing of ARISA and TRFLP data. It is designed to eliminate the laborious task of manually classifying community fingerprints in microbial community studies.}
\usage{RiboSort(data, dataformat = "standard", dye="B", output = "proportions", zerorows = FALSE, repeats = 1, mergerepeats = "none")}
\arguments{
  \item{data}{A character vector listing the filenames of each profile to be sorted. Each file listed must be stored in the current working directory to be detectable by R.}
  \item{dataformat}{The format of profiles listed contained in \code{data}. This must be one of "standard", "beckman", "abisingle" or "abimultiple".}
  \item{dye}{An indicator code identifying the dye that corresponds to the primer used. The argument \code{dye} only comes into effect when \code{dataformat} is "beckman" or "abisingle". The Beckman sequencer files generally have a dye indicator of the form "D1", "D4", etc., whereas the ABIsingle files generally just have a single letter indicator, eg. "B", "R", "Y", etc.}
  \item{output}{Specifies whether the object produced by RiboSort will contain abundances (as given in the original profiles) or the relative proportions of abundances in each profile. \code{output} must be one of "abundances" or "proportions". Statistical analysis is usually performed on the relative proportions of abundance, thus the default is "proportions".}
  \item{zerorows}{A logical indicating whether zerorows are to be kept in the output or not. A zerorow refers to a ribotype not detected in any of the profiles submitted. When \code{FALSE}, the default, zerorows are deleted from the output. When \code{TRUE}, the zerorows remain in the output.}
  \item{repeats}{The number of repeat profiles taken from each sample. If \code{repeats} is greater than 1, profiles listed in \code{data} must be in an order, such that all repeats from a particular sample are listed adjacent to one another. For example, if there were two repeats and five samples, Samp1repeat1 would be listed first followed by Samp1repeat2, Samp2repeat1, Samp2repeat2, etc.}
  \item{mergerepeats}{The method of merging a number of repeat profiles from the same sample into a single composite profile for that sample. \code{mergerepeats} must be one of "none", "presentinall", "presentintwo" or "presentinone".}
}
\details{
The number of repeats must be constant for all samples. If this is not the case, submit \code{repeats = 1} to obtain a sorted output of all profiles, and proceed to manually merge repeat profiles.

The \code{mergerepeats="none"} option indicates that profiles are not to be merged. To merge repeat profiles taken from the same sample into a single composite profile, there are three methods. The first of these, "presentinall", specifies that the composite profile only contains ribotypes detected in all of the repeat profiles. Thus, ribotypes present in less than all of the repeat profiles, are not included in the final composite profile. 

The second method, "presentintwo", specifies that the composite profile only contain ribotypes detected in at least two of the repeat profiles. Finally, the "presentinone" method indicates that all ribotypes detected, even those only present in one repeat profile, are included in the composite profile. This is the default option.

Following the execution of any of the three merging methods, the composite profile produced for a sample will be named according to the first repeat profile listed in \code{data} for that sample.

Composite profile abundances are determined by averaging the relative proportions present in repeat profiles.

When \code{repeats} is greater than 1 and  \code{mergerepeats="none"}, no profile merging will occur and all profiles submitted will be present in the output.
}
\value{A matrix with columns representing individual sample profiles, and assigned and aligned ribotypes in its rows.}
\author{Una Scallan, \email{una.scallan@ucd.ie}}
\section{Warning }{If an error appears stating 'subscript out of bounds', the most likely cause of this error is the submission of an inappropriate dye indicator for the given dataformat. Check the arguments submitted, ensuring the \code{dye} specified is that associated with your samples.}
\seealso{\code{\link{samplesMDS}}, \code{\link{speciesMDS}}}
\examples{
#Note: The four examples below make use of several demonstration datafiles provided with the RiboSort package. The line beginning file.copy(...) in each example, simply copies the neccessary demonstration datafiles from the RiboSort package into your current working directory. 

#Example 1: Data in standard profile format. 6 profiles, 2 repeats per sample.
#We specify that the RiboSort object x is to contain proportions, and the 
#present-in-one method is used to merge repeats.

data(demo.standard)
demo.standard
file.copy(system.file(demo.standard, package = "RiboSort"), ".")
x=RiboSort(data=demo.standard, dataformat="standard", output="proportions", 
	     zerorows=FALSE, repeats=2, mergerepeats="presentinone")
x[1:10,]		#prints the first 10 rows of the matrix x. 
samplesMDS(x)		#produces an MDS plot of the three samples.

#----------------------------------------------------------------------------------
#Example 2: Data in Beckman Sequencer format. 4 profiles, 1 repeat per sample.
#We specify that the RiboSort object y is to contain abundances.

data(demo.beckman)
demo.beckman
file.copy(system.file(demo.beckman, package = "RiboSort"), ".")
y=RiboSort(data=demo.beckman, dataformat="beckman", dye="D4", output="abundances",
           zerorows=FALSE, repeats=1, mergerepeats="none")
y[1:10,]		#prints the first 10 rows of the matrix y. 
samplesMDS(y)		#produces an MDS plot of the four samples.


#----------------------------------------------------------------------------------
#Example 3 : Data from ABI Sequencer files, each containing a single profile. 
#5 profiles, 1 repeat per sample.
#We specify that the RiboSort object z will contain abundances.

data(demo.abi.single)
demo.abi.single
file.copy(system.file(demo.abi.single, package = "RiboSort"), ".")
z=RiboSort(data=demo.abi.single, dataformat="abisingle", dye="B", output="abundances",
           zerorows=FALSE, repeats=1, mergerepeats="none")
z[1:10,]		#prints the first 10 rows of the matrix z. 
samplesMDS(z)		#produces an MDS plot of the five samples.



#----------------------------------------------------------------------------------
#Example 4: Data from an ABI Sequencer file that contains multiple profiles. 
#8 profiles, 1 repeat per sample.
#We specify that the RiboSort object q will contain abundances.

data(demo.abi.multiple)
demo.abi.multiple
file.copy(system.file(demo.abi.multiple, package = "RiboSort"), ".")
q=RiboSort(data=demo.abi.multiple, dataformat="abimultiple", output="abundances",
           zerorows=FALSE, repeats=1, mergerepeats="none")
q[1:10,]		#prints the first 10 rows of the matrix q. 
samplesMDS(q)		#produces an MDS plot of the eight samples.

}

\keyword{manip}
\keyword{multivariate}


