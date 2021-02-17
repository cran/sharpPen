testfun <-
function (x, k) 
{
    if (missing(k)) {
        6 * x + 3 * sin(x * (4 * pi)) + 5 * cos(x * (4 * pi))
    }
    else {
        2 - 5 * x + 2.5 * k * exp(-200 * k * (x - 0.5)^2)
    }
}
