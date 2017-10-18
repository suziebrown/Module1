# Define the class "lars"
# a list containing coefficients at each step, residuals at each step, L1 norm of coefficients at each step, coordinate updated at each step, method used

setClass("lars", slots=list(beta="matrix",resid="matrix",t="numeric" ,j="integer",method="character"))
