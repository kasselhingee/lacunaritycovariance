\name{balcattapark_coarse}
\docType{data}
\alias{balcattapark_coarse}
\title{Vegetation Mask of Balcatta Park}
\description{
 A vegetation map of Balcatta Park, Perth extracted from the 2007 or 2009** Urban Monitor data (see paper referenced below and report **). 
This data has been aggregated to a lower resolution of 0.8m per pixel from the original 0.2m (using the \code{raster} package's \code{aggregate} and a maximum function).
 
 A polygon approximating the boundary of the park is also included.

 The coordinate reference system for both vegetation map and polygon is in GDA94 / MGA Zone 50.
}
\usage{balcattapark_coarse}
\format{The vegetation map is stored in a binary mask in spatstat's owin format.
 The boundary polygon is also an owin object.}
\source{CSIRO Urban Monitor Project}
\references{
[1] P. Caccetta, S. Collings, A. Devereux, K. Hingee, D. McFarlane, A. Traylen, X. Wu, and Z.-S. Zhou, “Monitoring land surface and cover in urban and peri-urban environments using digital aerial photography,” International Journal of Digital Earth, vol. 0, no. 0, pp. 1–19, Jul. 2015.

[2] Caccetta, P., Collings, S., Devereux, A., Hingee, K., McFarlane, D., Traylen, A., Wu, X., and Zhou, Z., Urban Monitor: Enabling effective monitoring and management of urban and coastal environments using digital aerial photography; Final Report – Transformation of aerial photography into digital raster information product. CSIRO. 2012. https://www.environment.gov.au/resource/urban-monitor-enabling-effective-monitoring-and-management-urban-and-coastal-environments
}
\keyword{datasets}
