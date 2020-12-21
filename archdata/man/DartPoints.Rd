\name{DartPoints}
\alias{DartPoints}
\docType{data}
\title{
Five dart point types from Fort Hood, Texas, U.S.A.
}
\description{
Metric and categorical measurements on 91 Archaic dart points recovered during surface surveys at Fort Hood, Texas representing five types.
}
\usage{data(DartPoints)}
\format{
  A data frame with 91 observations on the following 17 variables.
  \describe{
    \item{\code{Name}}{Dart point type: \code{Darl}, \code{Ensor}, \code{Pedernales}, \code{Travis}, \code{Wells}}
    \item{\code{Catalog}}{Fort Hood catalog number}
    \item{\code{TARL}}{Texas Archeological Research Laboratory site number}
    \item{\code{Quad}}{Fort Hood Quad}
    \item{\code{Length}}{Maximum Length (mm)}
    \item{\code{Width}}{Maximum Width (mm)}
    \item{\code{Thickness}}{Maxmimum Thickness (mm)}
    \item{\code{B.Width}}{Basal width (mm)}
    \item{\code{J.Width}}{Juncture width (mm)}
    \item{\code{H.Length}}{Haft element length (mm)}
    \item{\code{Weight}}{Weight (gm)}
    \item{\code{Blade.Sh}}{Blade shape: \code{E - Excurvate}, \code{I - Incurvate}, \code{R - Recurvate}, \code{S - Straight}}
    \item{\code{Base.Sh}}{Base shape: \code{E - Excurvate}, \code{I - Incurvate}, \code{R - Recurvate}, \code{S - Straight}}
    \item{\code{Should.Sh}}{Shoulder shape: \code{E - Excurvate}, \code{I - Incurvate}, \code{S - Straight}, \code{X - None}}
    \item{\code{Should.Or}}{Shoulder orientation: \code{B - Barbed}, \code{H - Horizontal}, \code{T - Tapered}, \code{X - None}}
    \item{\code{Haft.Sh}}{Shape lateral haft element \code{A - Angular}, \code{E - Excurvate}, \code{I - Incurvate}, \code{R - Recurvate}, \code{S - Straight}}
    \item{\code{Haft.Or}}{Orientation lateral haft element: \code{C - Concave}, \code{E - Expanding}, \code{P - Parallel}, \code{T - Contracting}, \code{V - Convex}}
  }
}
\details{
Measurements on five types of dart points from Fort Hood in central Texas (Darl, Ensor, Pedernales, Travis, and Wells). The points were recovered during 10 different pedestrian survey projects during the 1980's and were classified and measured by H. Blaine Ensor using the system created by Futato (1983) as described in Carlson, S., et al 1987, pp 51-70 and Appendices 4 and 7.
}
\source{
Fort Hood Projectile Points. Electronic database compiling the results of multiple surface surveys at Fort Hood in the possession of David L. Carlson, Department of Anthropology, Texas A&M University, College Station, TX. The artifacts are curated at Fort Hood, TX by the Cultural Resources Branch of the Directorate of Public Works.
}
\references{
Carlson, David L. 2017. \emph{Quantitative Methods in Archaeology Using R}. Cambridge University Press, pp 52-60, 99-103, 106-107, 109-115, 148-157, 182-185, 198-211.

Carlson, S. B., H. B. Ensor, D. L. Carlson, E. A. Miller, and D E. Young. 1987. Archaeological Survey at Fort Hood, Texas Fiscal Year 1984. \emph{United States Army Fort Hood. Archaeological Resource Management Series, Research Report Number 14}.

Futato, E. M. 1983. Projectile Point Morphology: Steps Toward a Formal Account. in Proceedings of the Thirty-fourth Southeastern Archaeological Conference, Lafayette, Louisiana, October 27-19, 1977. \emph{Southeastern Archaeological Conference. Bulletin} 21: 38--81.

}
\examples{
data(DartPoints)
boxplot(Length~Name, DartPoints)
plot(Width~Length, DartPoints, pch=as.numeric(Name), main="FOrt Hood Dart Points")
legend("topleft", levels(DartPoints$Name), pch=1:5)
}
\keyword{datasets}
