#CP 03/10/2017 I want to keep the quantiles from wmr.boot, so this is just a slight improvement from the first mvcwt.wmr.boot function (Keitt 2014)

wmr.boot_perso=function (w, smoothing = 1, reps = 1000, mr.func = "wmr",aquantile=c(0.025,0.975)) 
{
#    require('foreach')
    require(mvcwt)
    mr.func = match.fun(mr.func)
    mr.obs = mr.func(w, smoothing = smoothing)
    with(w, {
        nloc = length(x)
        nvars = dim(z)[3]
        nscales = length(y)
        exports = c("reps", "wmr", "lagMat", "regularize", "mr.func", 
            "Gauss", "mr.obs", "nscales", "smoothing")
        flibs = c("mvcwt")
        mr.obs$z.boot = foreach(i = 1:nscales, .combine = c, 
            .export = exports, .packages = flibs) %dopar% {
            mr.boot = foreach(j = 1:reps, .combine = cbind, .inorder = FALSE) %dopar% 
                {
                  rphase = t(array(runif(nvars, -pi, pi), dim = c(nvars, 
                    nloc)))
                  zp = z[, i, , drop = FALSE] * complex(argument = rphase)
                  as.vector(mr.func(list(x = x, y = y[i], z = zp), 
                    smoothing = smoothing)$z)
                }
            res = foreach(j = 1:nloc, .combine = c) %dopar% {
                quantile(ecdf(mr.boot[j, ]),aquantile)
            }
            res
        }
        dim(mr.obs$z.boot) = c(length(aquantile),length(x), length(y))
        return(mr.obs)
    })
}

