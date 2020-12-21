
#Fucnion para simular patron individual de una especie focal en un patron multivariado, dejando el resto de especies en sus coordenadas
# pp: patron multivariado
# mimark: nombre de la especie (marca) que se quiere simular
# simulate: lista de patrones individuales simulados o expresion para simular el patron individual
#nsim: numero de simulaciones
# nmin: numero minimo de individuos esperados en el patron simulado.

multifocalsimulator <- function(pp, mimark, simulate,nsim=99,nmin=NULL){
    lista.sim<-list()
    pp.sp<-split(pp)
    if(is.list(simulate)){
       for (i in 1:nsim){
       progressreport(i,nsim)
         pp.sp[[mimark]]<-simulate[[i]] 
         split(pp)<-pp.sp
	 lista.sim[[i]]<- pp
	 }
    }
    if(is.expression(simulate)){
       for (i in 1:nsim){
       progressreport(i,nsim)
         pp.sp[[mimark]]<-eval(simulate) 
	 if(!is.null(nmin)){ # control para que los simulados no tengan menos de "nmin" puntos
	    while(pp.sp[[mimark]]$n<nmin){pp.sp[[mimark]]<-eval(simulate) }
	    }
         split(pp)<-pp.sp
	 lista.sim[[i]]<- pp
	 }
    }
    
    return(lista.sim)
}