
test_that("gower.agree(), confint(), and influence() work",
{
    data = matrix(c(1,2,3,3,2,1,4,1,2,NA,NA,NA,
                    1,2,3,3,2,2,4,1,2,5,NA,3,
                    NA,3,3,3,2,3,4,2,2,5,1,NA,
                    1,2,3,3,2,4,4,1,2,5,1,NA), 12, 4)
    set.seed(12)
    fit = gower.agree(data)
    mu.hat = round(fit$mu.hat, 3)
    names(mu.hat) = NULL
    expect_equal(mu.hat, 0.818)
    ci = round(confint(fit), 3)
    names(ci) = NULL
    expect_equal(ci[1], 0.555)
    expect_equal(ci[2], 0.971)
    ci = round(confint(fit, level = 0.99, type = "HPD"), 3)
    names(ci) = NULL
    expect_equal(ci[1], 0.554)
    expect_equal(ci[2], 0.99)
    inf = influence(fit, units = 6)
    dfbeta = round(inf$dfbeta.units, 3)
    names(dfbeta) = NULL
    expect_equal(dfbeta, -0.082)
})

