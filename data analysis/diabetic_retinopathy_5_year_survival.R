library(plotly)

# Load diabetic retinopathy study in wide format
read.csv("retinopathy_data.csv")

# Compute the pseudo-observations for the 5-year (60 month) survival probability based on the Dabrowska estimator
t0 <- c(60,60)
obs <- PO_func_dabrowska(obs,t0)

#fit a regression model using GEE with a logit link
fit <- geese(PO~age+mean_risk+type,
             id=id, data=obs,scale.fix=TRUE,family=gaussian,
             jack=TRUE, mean.link="logit",corstr="independence")
betas <- fit$beta
summary(fit)

des_mat <- model.matrix(PO~age+mean_risk+type, data = obs)
pred_values <- data.frame('Estimated survival'=exp(des_mat%*%betas)/(1+exp(des_mat%*%betas)),type=obs$type, Z=obs$age, X=obs$mean_risk)
pred_values <- pred_values[order(pred_values$Z,pred_values$X),]
p <- plot_ly(x=pred_values$Z, y=pred_values$X, z=pred_values$Estimated.survival, type="scatter3d", mode="markers",color = pred_values$type)
p <- layout(p, scene = list(xaxis = list(title = "Age at diagnosis", range = c(0,60)), yaxis = list(title = "Risk score", range = c(6,12)), zaxis = list(title = "5 year survival",range = c(0,0.55))))
print(p)
