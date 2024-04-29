install.packages("data.table")
library("data.table")

dta = "http://theta.edu.pl/wp-content/uploads/2018/03/dane_wino.txt"

dta <- readline(prompt = "Podaj ścieżkę do pliku: ")
data = read.table(dta, header = FALSE, sep = ',', dec = '.')
data = data[ , -1]
data <- na.omit(data)
head(data)


repeat{
    hipoteza <- as.numeric(readline(prompt = "Wybierz hipotezę: \n 1. Two.sided \n 2. Greater \n 3. Less \n Wybrana hipoteza:"))
        if(hipoteza == 1 | hipoteza == 2 | hipoteza == 3){
            break
    }
    else{
        print("Wybrałeś złą hipotezę.")
    }
}


#Test t-studenta dla jednej próby
nasz_t_test = function(proba, mu0, alternatywa) {
    mean_x = mean(proba)
    sd_x = sd(proba)
    n_x = length(proba)
    df_x = n_x - 1
    t = sqrt(n_x)*(mean_x-mu0)/sd_x
    print()
    if (alternatywa == 1) {
        p_wartosc = pt(t, df_x)
    } else if (alternatywa == 2) {
        p_wartosc = pt(t, df_x, lower.tail = FALSE)        
    } else {
        p_wartosc = 2*pt(abs(t), df_x, lower.tail = FALSE)
    }
    list(statystyka = t, p = p_wartosc)
}



    
#Test t=studenta dla dwóch prób niezależnych 
nasz_t_test2 = function(proba, proba2, alternatywa) {    
    mean_x = mean(proba)
    var_x = var(proba)
    n_x = length(proba)
    mean_y = mean(proba2)
    var_y = var(proba2)
    n_y = length(proba2)
    df = n_x + n_y - 2
    t = (mean_x-mean_y)/sqrt((n_x-1)*var_x+(n_y-1)*var_y)*sqrt((n_x*n_y*df) / (n_x+n_y))
    if (alternatywa == 1) {
        p_wartosc = pt(t, df)
    } else if (alternatywa == 2) {
        p_wartosc = pt(t, df, lower.tail = FALSE)        
    } else {
        p_wartosc = 2*pt(abs(t), df, lower.tail = FALSE) 
    }
    list(statystyka = t, p = p_wartosc)
}

#Test t-studenta dla dwóch prób zależnych
nasz_t_test3 = function(proba, proba2, alternatywa) {
    d = proba - proba2
    mean_d = mean(d)
    sd_d = sd(d)
    n_d = length(d)
    df_d = n_d - 1
    mu0 = 0
    t = (mean_d - mu0) / (sd_d / sqrt(n_d))
    if (alternatywa == 1) {
        p_wartosc = pt(t, df_d)
    } else if (alternatywa == 2) {
        p_wartosc = pt(t, df_d, lower.tail = FALSE)        
    } else {
        p_wartosc = 2*pt(abs(t), df_d, lower.tail = FALSE) 
    }
    list(statystyka = t, p = p_wartosc)
}


#Regresja liniowa
regresja_l <- function(y, int, ...) {
    lista_arg <- list(...)
     for (i in 1:length(lista_arg)) {
        print(class(lista_arg[[i]]))
    }
    
    if (int == 0) {
        X = as.matrix(do.call(cbind, lista_arg))
    } else {
        X = as.matrix(do.call(cbind, c(list(rep(1, length(y))), lista_arg)))
    }
    y = as.matrix(y)
    
    beta = solve(t(X)%*%X)%*%t(X)%*%y
    print(paste("number beta: ", length(beta)))
    if (int == 0) { 
        z <- beta[1] * lista_arg[[1]]
        for (i in 2:length(lista_arg)) {
          z <- z - beta[i] * lista_arg[[i]] }
        res <- y - z
      } else {
        # Komentarz: e = beta[n]*x[n-1]-beta[n-1]*x[n-2]-beta[n-2]*x[n-3]...
        e <- beta[2] * lista_arg[[1]]
        for (i in 2:length(lista_arg)) {
          e <- e - beta[i + 1] * lista_arg[[i]]
        }
        res <- y - beta[1] - e
      }

    S2e = t(res)%*%res / (nrow(y) - ncol(X))
    Se = sqrt(S2e)

    seBeta = sqrt(diag(c(S2e) * solve(t(X)%*%X)))
    t = beta / seBeta

    p.value.t = 2*pt(abs(t), nrow(y) - ncol(X), lower.tail = F)

    if (int == 0) {
        R2 = 1- t(res)%*%res / sum((y^2))
    } else {
        R2 = 1- t(res)%*%res / sum((y-mean(y))^2)
    }

    if (int == 0) {
        F = R2 / (1-R2) * (nrow(y) - ncol(X)) / (ncol(X))
        p.value.f = pf(F, ncol(X), nrow(y) - ncol(X), lower.tail = FALSE)
    } else {
        F = R2 / (1-R2) * (nrow(y) - ncol(X)) / (ncol(X) - 1)    
        p.value.f = pf(F, ncol(X) - 1, nrow(y) - ncol(X), lower.tail = FALSE)
    }
    
    results <- list(est = beta, p1 = p.value.t, r2 = R2, f = F, p2 = p.value.f)
    
}


#Analiza wariancji

data_anova = read.table(dta, header = FALSE, sep = ',', dec = '.')
data_anova <- na.omit(data_anova)

nasza_anova = function(data) {
  y = data[, 2:14] 
  alpha = data[, 1]  

  k = sort(unique(alpha))

  n.i = numeric(length(k))
  y.i = vector("list", length(k))  
  for (i in 1:length(k)) {
    n.i[i] = sum(alpha == k[i])
    y.i[[i]] = colMeans(y[alpha == k[i], ])
  }

  m.y = colMeans(y)

  SSb = sum(sapply(1:length(k), function(i) n.i[i] * sum((y.i[[i]] - m.y)^2)))
  S2b = SSb / (length(k) - 1)

  SSt = sum(colSums((sweep(y, 2, m.y, "-"))^2))
  S2t = SSt / (length(y) - 1)

  SSw = SSt - SSb
  S2w = SSw / (length(y) - length(k))

  F = S2b / S2w
  p = pf(F, length(k) - 1, length(y) - length(k), lower.tail = FALSE)

  #Wyszukiwanie grup o znacząco różniących się średnich
  diff_groups = c()
  for (i in 1:(length(k)-1)) {
    for (j in (i+1):length(k)) {
      if (max(abs(y.i[[i]] - y.i[[j]])) > sqrt(S2w / n.i[i] + S2w / n.i[j])) {
        diff_groups = c(diff_groups, paste(k[i], k[j], sep = "-"))
      }
    }
  }

  list(statistics = F, p.value = p, different_groups = diff_groups)
}



repeat{
    test <- readline(prompt = "Wybierz test: \n 1. Test t-studenta dla jednej próby \n 2. Test t-studenta dla dwóch prób niezależnych \n 3. Test t-studenta dla dwóch prób zależnych \n 4. Regresja liniowa \n 5. Analiza wariancji (ANOVA) \n Wybrany test: ")
    if(test == "1"){
        repeat{
            x <- readline(prompt = "Wybierz kolumnę:")
            x <- as.numeric(x)
            x <- as.numeric(data[ , x])
            srednia_x <- mean(x)
            print(srednia_x)
            break
        }
        nasz_t_test_result <- nasz_t_test(x, srednia_x, hipoteza)
        print(nasz_t_test_result)
        break
    }
    else if(test == "2"){
        repeat{
            x <- readline(prompt = "Wybierz kolumnę:")
            x <- as.numeric(x)
            x <- as.numeric(data[ , x])
            y <- readline(prompt = "Wybierz kolumnę2:")
            y <- as.numeric(y)
            y <- as.numeric(data[ , y])
            srednia_x <- mean(x)
            srednia_y <- mean(y)
            print(srednia_x)
            print(srednia_y)
            break
        }
        nasz_t_test2_result <- nasz_t_test2(x, y, hipoteza)
        print(nasz_t_test2_result)
        break
    } 
    else if(test == "3"){
        repeat{
            x <- readline(prompt = "Wybierz kolumnę:")
            x <- as.numeric(x)
            x <- as.numeric(data[ , x])
            y <- readline(prompt = "Wybierz kolumnę2:")
            y <- as.numeric(y)
            y <- as.numeric(data[ , y])
            srednia_x <- mean(x)
            srednia_y <- mean(y)
            print(srednia_x)
            print(srednia_y)
            break
        }
        nasz_t_test3_result <- nasz_t_test3(x, y, hipoteza)
        print(nasz_t_test3_result)
        break
    }
    else if(test == "4"){
        liczba_zmiennych = 0
        value1 <- ""

        # Sprawdzenie argumentu
        while (TRUE) {
          value1 <- readline(prompt = "Ile zmiennych niezależnych chcesz użyć w modelu:")
          if (is.na(as.numeric(value1))) {
            cat("Blad, podaj wartośc numeryczną: .\n")
          } else {
            break
          }
        }

        liczba_zmiennych <- as.numeric(value1)
        
        #Wybór zmiennych niezależnych
        niezalezne <- c()
        for (i in 1:liczba_zmiennych) {
            x <- readline(prompt = "Wybierz kolumnę:")
            x <- as.numeric(x)
            temp = as.numeric(data[ , x])
            niezalezne <- append(niezalezne,list(temp))
        }
        y <- readline(prompt = "Wybierz zmienną zależną (nr kolumny):")
        y <- as.numeric(y)
        y <- as.numeric(data[ , y])
        #Intercept 
        intercept <- readline(prompt = "Wybierz wartość dla argumentu intercept:")
        intercept <- as.numeric(intercept)
        result <- do.call(regresja_l, c(list(y, intercept), niezalezne))
        print(result)
        break
    }
    else if(test == "5"){
        result <- nasza_anova(data_anova)
        print(result)
        break
    }
    else{
        print("Wybrałeś zły test")
    }     
}




library(MASS)

mme_2 = function(y, X, Z, A, sigma_a, sigma_e) {
  alpha = sigma_e / sigma_a
  invA = ginv(A)
C = rbind(cbind(t(X) %*% X, t(X) %*% Z),
          cbind(t(Z) %*% X, t(Z) %*% Z + alpha * t(Z) %*% invA %*% Z))

  rhs = rbind(t(X) %*% y, t(Z) %*% y)
  invC = ginv(C)
  estimators = invC %*% rhs
  list(C = C, est = estimators)
}

EM_2 = function(y, X, Z, A, sigma_a, sigma_e) {
  n = nrow(X)
  p = ncol(X)
  q = nrow(A)
  
  t = 1 
  tmp = 0.1 
  
  while (tmp > 0.00001) {
    mme_new = mme_2(y, X, Z, A, sigma_a, sigma_e)
    C_new = ginv(mme_new$C)
    Ck = C_new[(p + 1):(p + q), (p + 1):(p + q)]
    mme2 = mme_new$est
    
    a = as.matrix(mme2[(p + 1):(p + q)])
    sigma_a_new = (t(a) %*% ginv(A) %*% a + sum(diag(ginv(A) %*% Ck)) * c(sigma_e)) / q
    
    res = as.matrix(y - X %*% as.matrix(mme2[1:p]) - Z %*% as.matrix(mme2[(p + 1):(p + q)]))
    X.tmp1 = cbind(X, Z) %*% C_new
    X.tmp2 = t(cbind(X, Z))
    sigma_e_new = (t(res) %*% res + sum(diag(X.tmp1 %*% X.tmp2)) * c(sigma_e)) / n
    
    tmp = max(abs(sigma_a - sigma_a_new), abs(sigma_e - sigma_e_new))
    sigma_a = sigma_a_new
    sigma_e = sigma_e_new
    
    t = t + 1
  }
  
  list(t = t, sigma_a = sigma_a, sigma_e = sigma_e)
}


y <- data_anova[, 1]  #pierwsza kolumna jako zmienna y
X <- as.matrix(data_anova[, -1])  # wszystkie kolumny oprócz pierwszej jako zmienne efektów stałych
Z <- matrix(1, nrow = nrow(data_anova), ncol = 3)  #zmienne efektów losowych
A <- diag(3)  #macierz dla efektów losowych
sigma_a <- 20
sigma_e <- 40

results <- EM_2(y, X, Z, A, sigma_a, sigma_e)

#estymatory efektów losowych i stałych
est_random <- results$est[(ncol(X) + 1):ncol(results$est)]
est_fixed <- results$est[1:ncol(X)]

#estymatory wariancji
est_sigma_a <- results$sigma_a
est_sigma_e <- results$sigma_e

print("Estimators of random effects:")
print(est_random)
print("Estimators of fixed effects:")
print(est_fixed)
print("Estimator of sigma_a:")
print(est_sigma_a)
print("Estimator of sigma_e:")
print(est_sigma_e)

wine_types <- as.factor(data_anova[, 1])

#ANOVA
anova_results <- summary(aov(data_anova[, 2] ~ wine_types, data = data_anova))

print(anova_results)

#znaczące różnice między grupami
if (anova_results$Pr[1] < 0.05) {
  print("There are significant differences between groups.")
  significant_groups <- paste(anova_results$rownames[1], "and", anova_results$rownames[-1][anova_results$Pr[-1] < 0.05])
  print(paste("Significant differences exist between:", significant_groups))
} else {
  print("There are no significant differences between groups.")
}



