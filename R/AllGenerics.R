setGeneric(
	name = "sim",
	def=function(object, params){standardGeneric("sim")}
)


setGeneric(
	name="setOntology<-", 
	def=function(object, value) {standardGeneric("setOntology<-")}
)


setGeneric(
	name="setOrganism<-", 
	def=function(object, value) {standardGeneric("setOrganism<-")}
)

setGeneric(
	name="setMethod<-", 
	def=function(object, value) {standardGeneric("setMethod<-")}
)

setGeneric(
	name="setCombineMethod<-", 
	def=function(object, value) {standardGeneric("setCombineMethod<-")}
)

setGeneric(
	name="computeIC", 
	def=function(params) {standardGeneric("computeIC")}
)


setGeneric(
	name="loadGOMap", 
	def=function(params) {standardGeneric("loadGOMap")}
)

setGeneric(
	name="loadICdata", 
	def=function(params) {standardGeneric("loadICdata")}
)

