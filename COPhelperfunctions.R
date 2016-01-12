theta.eps = 0.01

ml.sample = function(u1, copula.typeL = GUMBELHAC) {
  d = density.3d(u1, Ltype = copula.typeL)
  mylist = list(thetas = d$thetas, ml = d$loglik, str = d$str, start = 0, end = 0)
  mylist
}

sample.u = function() {
  u = rcopula.GumbelNested(samplesize, tau2theta(c(tau1, tau2)))
  a = u[(brake.point + 1):(samplesize), 1]
  u[(brake.point + 1):(samplesize), 1] = u[(brake.point + 1):(samplesize), 3]
  u[(brake.point + 1):(samplesize), 3] = a
  u
}

sample.u.changeInFirstParam = function() {
  u1 = rcopula.GumbelNested(brake.point, tau2theta(c(tau1_first, tau2_general)))
  u2 = rcopula.GumbelNested(samplesize - brake.point, tau2theta(c(tau1_second, tau2_general)))
  rbind(u1, u2)
}

sample.u.changeInSecondParam = function() {
  u1 = rcopula.GumbelNested(brake.point, tau2theta(c(tau1_general, tau2_first)))
  u2 = rcopula.GumbelNested(samplesize - brake.point, tau2theta(c(tau1_general, tau2_second)))
  rbind(u1, u2)
}

make.crit.val = function(type.of.change = 1, num.sample = 0, u = NA, copula_typeL = "gumbel") {
  Ltype = if (copula_typeL == "gumbel") {
    GUMBELHAC
  } else if (copula_typeL == "clayton") {
    CLAYTONHAC
  }
  if (sum(is.na(u)) != FALSE) {
    if (type.of.change == 1) {
      u = sample.u.changeInFirstParam()
    }
    if (type.of.change == 2) {
      u = sample.u.changeInSecondParam()
    }
    if (type.of.change == 3) {
      u = sample.u()
    }
  }
  A = rep(0, 5)
  T1 = skip  
  MLI.old = 0
  k = 1
  while (T1 < samplesize) {
    MLI = ml.sample(u[(T1 - trunc(m0LCP * cLCP^(k + 1))):T1, ], Ltype)
    MLI$start = (T1 - trunc(m0LCP * cLCP^(k + 1)))
    MLI$end = T1
    maxtstat = -1e+05
    for (Lt in (T1 - trunc(m0LCP * cLCP^k)):(T1 - trunc(m0LCP * cLCP^(k - 1) - 1))) {
      tstat = ml.sample(u[(T1 - trunc(m0LCP * cLCP^(k + 1))):Lt, ], Ltype)$ml + ml.sample(u[Lt:T1, ], 
              Ltype)$ml - MLI$ml
      if (maxtstat < tstat) {
        maxtstat = tstat
      }
    }
    a = trunc(sort(c(taus[which.min(abs(theta2tau(MLI$thetas, type = copula_typeL)[1] - taus))], 
                     taus[which.min(abs(theta2tau(MLI$thetas, type = copula_typeL)[2] - taus))])) * 10)
    if (maxtstat < qx[paste(a[1], a[2], sep = ""), k]) {
      if (k != Ak_max) {
        k = k + 1
        MLI.old = MLI
      } else {
        A = rbind(A, c((MLI.old$end - MLI.old$start), conv.str(MLI.old$str), MLI.old$theta[1], MLI.old$theta[2], 
                       MLI.old$ml/(MLI.old$end - MLI.old$start)))
        T1 = T1 + 1
        print(paste(" T = ", (T1 - skip), " L = ", (MLI.old$end - MLI.old$start), sep = ""))
        k = 1
      }
    } else {
      MLI = ml.sample(u[(T1 - trunc(m0LCP * cLCP^(k - 1))):T1, ])
      MLI$start = (T1 - trunc(m0LCP * cLCP^(k - 1)))
      MLI$end = T1
      MLI.old = MLI
      A = rbind(A, c((MLI.old$end - MLI.old$start), conv.str(MLI.old$str), MLI.old$theta[1], MLI.old$theta[2], 
                     MLI.old$ml/(MLI.old$end - MLI.old$start)))
      T1 = T1 + 1
      print(paste(" T = ", (T1 - skip), " L = ", (MLI.old$end - MLI.old$start), sep = ""))
      k = 1
    }
  }
  A = A[-1, ]
  colnames(A) = c("length", "structure", "theta1", "theta2", "ML")
  A
}

density.3d = function(u1, Ltype = GUMBELHAC, only.cdf = F) {
  if ((Ltype == GUMBELHAC) || (Ltype == GUMBELAC)) {
    gen_type = "gumbel"
  }
  if ((Ltype == CLAYTONHAC) || (Ltype == CLAYTONAC)) {
    gen_type = "clayton"
  }
  if ((Ltype == GUMBELHAC) || (Ltype == CLAYTONHAC)) {
    k1 = cor(u1, method = "kendall") - diag(rep(1, 3))
    k1 = k1[upper.tri(k1)]
    param = max(k1)
    max.k1.ind = which(k1 == param)
    if (length(max.k1.ind) > 1) {
      max.k1.ind = max.k1.ind[1]
    }
    if (max.k1.ind == 1) {
      param = c(param, 1, 2)
    } else if (max.k1.ind == 2) {
      param = c(param, 1, 3)
    } else if (max.k1.ind == 3) {
      param = c(param, 2, 3)
    }
    theta1 = max(min.cop(gen_type), tau2theta(param[1], gen_type))
    theta = c(max(min.cop(gen_type), tau2theta(cor(u1[, ret.rest(sum(param[2:3]))], 
            cop2d(u1[, param[2]], u1[, param[3]], theta1, gen_type), method = "kendall"), gen_type)), theta1)
    if (theta[2] < theta[1]) {
      theta[1] = theta[2] - theta.eps/2
    }
    if (abs(theta[1] - theta[2]) < theta.eps) {
      if (only.cdf == F) {
        dens = cop.123.density(u1[, param[2]], u1[, param[3]], u1[, ret.rest(sum(param[2:3]))], theta[1], gen_type)
      } else {
        dens = NA
      }
      cdf.cop = cop.123.cdf(u1[, param[2]], u1[, param[3]], u1[, ret.rest(sum(param[2:3]))], theta[1], gen_type)
      
    } else {
      if (only.cdf == F) {
        dens = cop.12.3.density(u1[, param[2]], u1[, param[3]], u1[, ret.rest(sum(param[2:3]))], theta[1], theta[2], 
                                gen_type)
      } else {
        dens = NA
      }
      cdf.cop = cop.12.3.cdf(u1[, param[2]], u1[, param[3]], u1[, ret.rest(sum(param[2:3]))], theta[1], theta[2], 
                             gen_type)
    }
    if (only.cdf == F) {
      Lloglik = sum(log(dens))
      loglik = Lloglik
    } else {
      Lloglik = NA
    }
    theta.matr = matrix(0, nrow = 3, ncol = 3)
    theta.matr[param[2], param[3]] = theta[2]
    theta.matr[param[3], param[2]] = theta[2]
    theta.matr[param[2], ret.rest(sum(param[2:3]))] = theta[1]
    theta.matr[ret.rest(sum(param[2:3])), param[2]] = theta[1]
    theta.matr[param[3], ret.rest(sum(param[2:3]))] = theta[1]
    theta.matr[ret.rest(sum(param[2:3])), param[3]] = theta[1]
    
    mylist = list(thetas = theta, density = dens, cdf = cdf.cop, 
                  str = paste("((", round(param[2]), ".", round(param[3]), ").", 
                              ret.rest(sum(param[2:3])), ")", sep = ""), 
                  loglik = Lloglik, thetaM = theta.matr, BIC = (-2 * Lloglik + 2 * log(2)))
  }
  if (Ltype == GUMBELAC) {
    gumbel.cop = gumbelCopula(1.5, dim = 3)
    fit = fitCopula(gumbel.cop, u1, method = "ml")
    gumbel.cop = gumbelCopula(fit@estimate, dim = 3)
    
    if (only.cdf == F) {
      dens = dcopula(gumbel.cop, u1)
    } else {
      dens = NA
    }
    mylist = list(thetas = fit@estimate, density = dens, cdf = pcopula(gumbel.cop, u1), 
                  loglik = fit@loglik, BIC = (-2 * fit@loglik + 2 * log(1)))
  }
  if (Ltype == CLAYTONAC) {
    clayton.cop = claytonCopula(1.5, dim = 3)
    fit = fitCopula(clayton.cop, u1, method = "ml")
    clayton.cop = claytonCopula(fit@estimate, dim = 3)
    
    if (only.cdf == F) {
      dens = dcopula(clayton.cop, u1)
    } else {
      dens = NA
    }
    mylist = list(thetas = fit@estimate, density = dens, cdf = pcopula(clayton.cop, u1), 
                  loglik = fit@loglik, BIC = (-2 * fit@loglik + 2 * log(1)))
  }
  if (Ltype == GAUSS) {
    cor.data = cor(qnorm(u1))
    cor.param = cor.data[lower.tri(cor.data)]
    normal.cop = normalCopula(cor.param, dim = 3, dispstr = "un")
    
    if (only.cdf == F) {
      dens = dcopula(normal.cop, u1)
      Lloglik = sum(log(dens))
    } else {
      dens = NA
      Lloglik = NA
    }
    
    mylist = list(thetas = cor.param, density = dens, cdf = pcopula(normal.cop, u1), 
                  loglik = Lloglik, BIC = (-2 * Lloglik + 2 * log(3)))
  }
  mylist
}

min.cop = function(Ltype) {
  if (Ltype == "gumbel") {
    1.01
  } else if (Ltype == "clayton") {
    0.01
  }
}

tau2theta = function(tau, type = "gumbel") {
  if (type == "gumbel") {
    1/(1 - tau)
  } else if (type == "clayton") {
    2 * tau/(1 - tau)
  } else if (type == "gauss") {
    sin(tau * pi/2)
  }
}

cop2d = function(x1, x2, Ltheta, Ltype = "gumbel") {
  phi(phi_m1(x1, Ltheta, Ltype) + phi_m1(x2, Ltheta, Ltype), Ltheta, Ltype)
}

phi = function(x, theta, type = "gumbel") {
  if (type == "gumbel") {
    exp(-x^{
      1/theta
    })
  } else {
    (theta * x + 1)^(-1/theta)
  }
}

phi_m1 = function(x, theta, type = "gumbel") {
  if (type == "gumbel") {
    (-log(x))^theta
  } else {
    (x^(-theta) - 1)/theta
  }
}

ret.rest = function(d) {
  if (d == 3) {
    3
  } else if (d == 4) {
    2
  } else if (d == 5) {
    1
  }
}

cop.12.3.density = function(u1, u2, u3, theta1, theta2, type = "gumbel") {
  if (type == "gumbel") {
    gumb.12.3.density(u1, u2, u3, theta1, theta2)
  } else if (type == "clayton") {
    clayton.12.3.density(u1, u2, u3, theta1, theta2)
  }
}

gumb.12.3.density = function(u1, u2, u3, theta1, theta2) {
  l1 = -log(u1)
  l2 = -log(u2)
  l3 = -log(u3)
  theta1m1 = (theta1 - 1)
  theta2m1 = (theta2 - 1)
  onedtheta1 = 1/theta1
  l12.2 = (l1)^theta2 + (l2)^theta2
  l12.2.t2 = l12.2^(1/theta2)
  c12 = (l12.2.t2)^theta1
  c12.3 = c12 + l3^theta1
  tc12.3 = (theta1m1 + c12.3^onedtheta1)
  
  -((l12.2.t2^(theta1 - 2) * (l1^theta2m1) * l12.2^(-2 + 1/theta2) * l2^theta2m1 * 
   (theta2m1 * (-l12.2.t2) * tc12.3 * c12.3 + l12.2.t2 * 
   (-(c12 * (theta1m1 * theta1 + 2 * theta1m1 * c12.3^onedtheta1 + c12.3^(2 * onedtheta1))) + 
   theta1m1 * tc12.3 * l3^theta1)) * c12.3^(onedtheta1 - 3) * l3^theta1m1)/
   (exp(c12.3^onedtheta1) * u1 * u2 * u3))
}

cop.12.3.cdf = function(x1, x2, x3, Ltheta1, Ltheta2, Ltype = "gumbel") {
  cop2d(cop2d(x1, x2, Ltheta2, Ltype), x3, Ltheta1, Ltype)
}

theta2tau = function(theta, type = "gumbel") {
  if (type == "gumbel") {
    1 - 1/theta
  } else if (type == "clayton") {
    theta/(2 + theta)
  } else if (type == "gauss") {
    2/pi * asin(theta)
  }
}

cop.123.density = function(u1, u2, u3, theta, type = "gumbel") {
  if (type == "gumbel") {
    gumb.123.density(u1, u2, u3, theta)
  } else if (type == "clayton") {
    clayton.123.density(u1, u2, u3, theta)
  }
}

gumb.123.density = function(u1, u2, u3, theta) {
  lu1 = -log(u1)
  lu2 = -log(u2)
  lu3 = -log(u3)
  lu1t = lu1^theta
  lu2t = lu2^theta
  lu3t = lu3^theta
  lu123t = lu1t + lu2t + lu3t
  
  (lu1^(-1 + theta) * lu2^(-1 + theta) * (1 - 3 * theta + 2 * theta^2 + 3 * 
  (-1 + theta) * (lu123t)^(1/theta) + lu123t^(2/theta)) * 
  lu123t^(-3 + 1/theta) * lu3^(-1 + theta))/(exp(lu123t^(1/theta)) * u1 * u2 * u3)
}

cop.123.cdf = function(x1, x2, x3, theta, type = "gumbel") {
  phi(phi_m1(x1, theta, type) + phi_m1(x2, theta, type) + phi_m1(x3, theta, type), theta, type)
}

clayton.12.3.density = function(u1, u2, u3, theta2, theta1) {
  u1t1u2t1 = -1 + u1^(-theta1) + u2^(-theta1)
  it1 = 1/theta1
  it2 = 1/theta2
  u1mt1 = u1^(-1 - theta1)
  u2mt1 = u2^(-1 - theta1)
  u3mt2 = u3^(-1 - theta2)
  u1t1u2t1it1 = u1t1u2t1^(-it1)
  shorter = -1 + u1t1u2t1it1^(-theta2) + u3^(-theta2)
  long = shorter^(-2 - it2)
  
  (-2 - it2) * (-1 - it2) * (theta2^2) * u1mt1 * u2mt1 * (u1t1u2t1^(-2 - 2 * it1)) * 
    (u1t1u2t1it1^(-2 - 2 * theta2)) * 
    u3mt2 * (shorter^(-3 - it2)) - (-1 - it2) * (-1 - theta2) * 
    theta2 * u1mt1 * u2mt1 * (u1t1u2t1^(-2 - 2/theta1)) * 
    (u1t1u2t1it1^(-2 - theta2)) * u3mt2 * long + (-1 - it1) * theta1 * (-1 - it2) * 
    theta2 * u1mt1 * u2mt1 * (u1t1u2t1^(-2 - it1)) * (u1t1u2t1it1^(-1 - theta2)) * u3mt2 * long
}

clayton.123.density = function(u1, u2, u3, theta) {
  (2 + 1 / theta) * theta * (1 + theta) * (u1^(-1 - theta)) * (u2^(-1 - theta)) *
    (u3^(-1 - theta)) * (-2 + u1^(-theta) + u2^(-theta) + u3^(-theta))^(-3 - 1 / theta)
}

emp.copula.2d = function(u1, u2, Lcdf.x) {
  n = length(u1)
  colSums(matrix(rowSums((matrix(rep(cbind(Lcdf.x[, 1], Lcdf.x[, 2]), each = n), ncol = 2) >= 
                           cbind(rep(Lcdf.x[, 1], n), rep(Lcdf.x[, 2], n)))), ncol = n) == 2)/n
}

emp.copula.self.2d = function(Lcdf.x) {
  emp.copula.2d(Lcdf.x[, 1], Lcdf.x[, 2], Lcdf.x)
}

alpha = 0.95
z_alpha = qnorm((1 + alpha)/2) 