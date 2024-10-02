\name{glmer_cb2_short}
\alias{glmer_cb2_short}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Logistic Regression Mixed Model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Function designed to parallelize glmer using dplyr's group() and summarize()
piping for each position independently in a BSA experiment. Function takes in
parameters from tibble and outputs rows for each result of the glmer (one
per coefficient).
}
\usage{
glmer_cb2_short(..., W, formula, outputlength = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
Define all parameters to used IN FORMULA here; do not include parameters which
are NOT included in the formula specifically. These will be called by name by
the formula within the function.
}
  \item{W}{
%%     ~~Describe \code{W} here~~
Name of column to be used as binomial weight to insert into glmer()
}
  \item{formula}{
%%     ~~Describe \code{formula} here~~
Text in quotes for input function (see glmer for details).
}
  \item{outputlength}{
%%     ~~Describe \code{outputlength} here~~
Number of rows expected for single output of each coefficient. If this is not
satisfied, that position will be skipped so that the pipe does not crash. Set to
4 as the default for formula "C ~ A + B + A*B + (1 | R)", as the intercept, A,
B, and interaction terms will each produce a coefficient.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
If using a simple linear model (no random effects), use glm_cb2_short() instead,
as it will call glm() instead of glmer().
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Column of z-scores for each coefficient in formula, in the order that is
returned by the glmer() function. Test a single position first without piping
to determine what labels should be used (see second example for adding labels
to summary).
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
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
# Make data
require(dplyr)
sampledf <- data.frame(Allele = c(1,2),
                       POS = 1,
                       Bulk = c(1,1,1,1,2,2,2,2),
                       Rep = c(1,1,2,2,1,1,2,2),
                       Count = rnorm(n = 8)
                       ) %>% mutate_if(is.numeric, as.factor) %>%
    mutate(Count = as.numeric(Count))

sampledf %>%
    group_by(POS) %>%
    summarize(zscore = glmer_cb2_short(Allele = Allele, Bulk = Bulk, Rep = Rep,
                                        W = Count,
                                        formula = "Allele ~ Bulk + (1 | Rep)",
                                        output = 2)) -> sampleoutput

sampledf %>%
    group_by(POS) %>%
    summarize(zscore = glmer_cb2_short(Allele = Allele, Bulk = Bulk, Rep = Rep,
                                        W = Count,
                                        formula = "Allele ~ Bulk + (1 | Rep)",
                                        output = 2),
              label = c("Intercept", "Bulk")) -> sampleoutput

################################################################################

## The function is currently defined as
function (..., W, formula, numgroups = FALSE, outputlength = 4,
    return = c("Z"))
{
    data <- list(...)
    require(lme4)
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
    glm_fit <- glmer(glm_formula, data = as.data.frame(data),
        weights = W, family = binomial)
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
