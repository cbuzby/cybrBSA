\name{glm_cb2_short}
\alias{glm_cb2_short}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Logistic Regression Model
}
\description{
Function designed to parallelize glmer using dplyr's group() and summarize() piping for each position independently in a BSA experiment. Function takes in parameters from tibble and outputs rows for each result of the glmer (one per coefficient).


}
\usage{
glm_cb2_short(..., W, formula, numgroups = FALSE, outputlength = 4, return = c("Z"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
Define all parameters to used IN FORMULA here; do not include parameters which are NOT included in the formula specifically. These will be called by name by the formula within the function.
}
  \item{W}{
Name of column to be used as binomial weight to insert into glm()
}
  \item{formula}{
Text in quotes for input function (see glm for details).
}
  \item{numgroups}{
Deprecated, set to false
}
  \item{outputlength}{
Number of rows expected for single output of each coefficient. If this is not satisfied, that position will be skipped so that the pipe does not crash. Set to 4 as the default for formula "C ~ A + B + A*B", as the intercept, A, B, and interaction terms will each produce a coefficient.
}
  \item{return}{
Deprecated, set to "Z" to return z-scores, but has no alternative
}
}
\details{
If using a mixed linear model (including random effects), use glmer_cb2_short() instead, as it will call glmer() instead of glm().
}
\value{
Column of z-scores for each coefficient in formula, in the order that is returned by the glmer() function. Test a single position first without piping to determine what labels should be used (see second example for adding labels to summary).
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Cassandra Buzby
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or standard data sets, see data().

## The function is currently defined as
function (..., W, formula, numgroups = FALSE, outputlength = 4,
    return = c("Z"))
{
    data <- list(...)
    if (is.null(W) || is.null(formula)) {
        stop("Weights (W) and formula must be provided")
    }
    glm_formula <- as.formula(formula)
    if (!all(names(data) \%in\% all.vars(glm_formula))) {
        stop("One or more variables in the formula are not provided as arguments")
    }
    for (i in all.vars(glm_formula)) {
        if (length(unique(as.data.frame(data)[, i])) < 2) {
            output <- rep(NA, outputlength)
            return(output)
        }
    }
    glm_fit <- glm(glm_formula, data = as.data.frame(data), weights = W,
        family = binomial)
    if (return \%in\% "Z") {
        output <- summary(glm_fit)$coefficients[((length(summary(glm_fit)$coefficients) *
            0.5) + 1):((length(summary(glm_fit)$coefficients) *
            0.75))]
    }
    if (length(output) == outputlength) {
        return(output)
    }
    else {
        return(rep(NA, outputlength))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
