\name{yeast}
\alias{yeast}

\title{Yeast Cell Cycle Data Set}
\description{ This is the Yeast Cell Cycle dataset that was originally analyzed in Chun and Keles (2010). It is also available in the package, \code{spls}.
 }
\usage{data(yeast)}
\format{
A list of two components.
\describe{
    \item{\code{x}:}{ ChIP-chip data. A matrix with 542 rows and 106 columns.}
    \item{\code{y}:}{ Cell cycle gene expression data. A matrix with 542 rows and 18 columns.}
  }
}

\details{ Matrix \eqn{y} is cell cycle gene expression data (Spellman et al., 1998) of 542 genes from an \eqn{\alpha} factor based experiment. Each column corresponds to mRNA levels measured at every 7 minutes during 119 minutes (a total of 18 measurements). Matrix \eqn{x} is the chromatin immunoprecipitation on chip (ChIP-chip) data of Lee et al. (2002) and it contains the binding information for 106 transcription
factors. See Chun and Keles (2010) for more details.
}
\source{ 
Lee, T.I., Rinaldi, N.J., Robert, F., Odom, D.T., Bar-Joseph, Z., Gerber, G.K., Hannett, N.M., Harbison, C.T., Thomson, C.M., Simon, I., Zeitlinger, J., Jennings, E.G., Murray, H.L., Gordon, D.B., Ren, B., Wyrick, J.J., Tagne, J.B., Volkert, T.L, Fraenkel, E, Gifford, D.K., and Young, R.A (2002). Transcriptional Regulatory Networks in Saccharomyces Cerevisiae. \emph{Science}, \bold{298}: 799-804.

Spellman, P.T., Sherlock, G., Zhang, M.Q., Iyer, V.R., Anders, K., Eisen, M.B., Brown, P.O., Botstein, D., and Futcher, B. (1998). Comprehensive Identification of Cell Cycle-Regulated Genes of the Yeast
Saccharomyces Cerevisiae by Microarray Hydrization. \emph{Molecular Biology of the Cell}, \bold{9}: 3273-3279.
}
\references{
Chun, H. and Keles, S. (2010). Sparse Partial Least Squares for Simultaneous Dimension Reduction and Variable Selection. \emph{Journal of the Royal Statistical Society - Series B}, \bold{72}: 3-25.
}