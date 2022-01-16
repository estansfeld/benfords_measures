# some background
# http://statistics.gov.scot/data/deaths-involving-coronavirus-covid-19
# https://medium.swirrl.com/using-r-to-analyse-linked-data-7225eefe2eb8
#
# training resources
# https://guides.statistics.gov.scot/article/22-querying-data-with-sparql
# https://guides.statistics.gov.scot/article/34-understanding-the-data-structure

library(SPARQL)
get_SG_CovidDeaths_Region_SparQL <- function() {
  
  # first we need to connect to the SPARQL end-point, which is
  endpoint <- 'http://statistics.gov.scot/sparql'
  
  # The SG data uses the schema and data dictionary defined below
  query <- 'PREFIX qb: <http://purl.org/linked-data/cube#>
    PREFIX data: <http://statistics.gov.scot/data/>
    PREFIX sdmxd: <http://purl.org/linked-data/sdmx/2009/dimension#>
    PREFIX defs: <http://statistics.gov.scot/def/dimension/>
    PREFIX mp: <http://statistics.gov.scot/def/measure-properties/>
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>

    SELECT ?periodname ?area ?areauri ?nDeaths 

    WHERE {
     ?obs qb:dataSet data:deaths-involving-coronavirus-covid-19 .
    	?obs sdmxd:refArea  ?areauri .
    	?obs sdmxd:refPeriod ?perioduri .
    	?obs mp:count ?nDeaths .
    
      ?areauri rdfs:label ?area .
    	?perioduri rdfs:label ?periodname .

}'
  return(SPARQL(endpoint,query))
}


get_SG_CovidDeaths_SparQL <- function() {
  
  # first we need to connect to the SPARQL end-point, which is
  endpoint <- 'http://statistics.gov.scot/sparql'
  
  # The SG data uses the schema and data dictionary defined below
  query <- 'PREFIX qb: <http://purl.org/linked-data/cube#>
    PREFIX data: <http://statistics.gov.scot/data/>
    PREFIX sdmxd: <http://purl.org/linked-data/sdmx/2009/dimension#>
    PREFIX defs: <http://statistics.gov.scot/def/dimension/>
    PREFIX mp: <http://statistics.gov.scot/def/measure-properties/>
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>

    SELECT ?periodname ?area ?areauri ?sex ?cause ?location ?age ?nDeaths 

    WHERE {
     ?obs qb:dataSet data:deaths-involving-coronavirus-covid-19 .
      ?obs defs:age ?ageuri .  
    	?obs defs:causeOfDeath ?causeuri .   
    	?obs defs:locationOfDeath ?locationuri .  
    	?obs sdmxd:refArea  ?areauri .
    	?obs sdmxd:refPeriod ?perioduri .
    	?obs defs:sex ?sexuri . 
    	?obs mp:count ?nDeaths .
    
      ?ageuri rdfs:label ?age .
    	?causeuri rdfs:label ?cause .   
    	?locationuri rdfs:label ?location .  
      ?areauri rdfs:label ?area .
    	?perioduri rdfs:label ?periodname .
    	?sexuri rdfs:label ?sex .

}'
  return(SPARQL(endpoint,query))
}


get_SG_CovidInfo_SparQL <- function() {
  
  # first we need to connect to the SPARQL end-point, which is
  endpoint <- 'http://statistics.gov.scot/sparql'
  
  # The SG data uses the schema and data dictionary defined below
  query <- 'PREFIX qb: <http://purl.org/linked-data/cube#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX sdmx: <http://purl.org/linked-data/sdmx/2009/concept#>
PREFIX data: <http://statistics.gov.scot/data/>
PREFIX sdmxd: <http://purl.org/linked-data/sdmx/2009/dimension#>
PREFIX mp: <http://statistics.gov.scot/def/measure-properties/>
PREFIX stat: <http://statistics.data.gov.uk/def/statistical-entity#>
SELECT * WHERE {
?indicator qb:dataSet data:coronavirus-covid-19-management-information ;
             sdmxd:refArea ?area ;
             sdmxd:refPeriod ?year ;
             mp:count ?ncount .
?year rdfs:label ?yearname .
  
?area stat:code ?areatype ;
      rdfs:label ?areaname .
?areatype rdfs:label ?areatypename .
}'
  return(SPARQL(endpoint,query))
}

get_SG_Population_SparQL <- function() {
  
  # first we need to connect to the SPARQL end-point, which is
  endpoint <- 'http://statistics.gov.scot/sparql'
  
  # The SG data uses the schema and data dictionary defined below
  query <- "PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
    PREFIX qb: <http://purl.org/linked-data/cube#>
    PREFIX sdmx: <http://purl.org/linked-data/sdmx/2009/dimension#>
    PREFIX dim: <http://statistics.gov.scot/def/dimension/>
    PREFIX measures: <http://statistics.gov.scot/def/measure-properties/>
    
SELECT ?periodname ?area ?age ?sex ?projection ?nPeople 

WHERE {
	?obs qb:dataSet <http://statistics.gov.scot/data/population-projections-2018-based> .
	?obs sdmx:refArea ?areauri .
	?obs sdmx:refPeriod ?perioduri .
	?obs dim:age ?ageuri .
  ?obs dim:sex ?sexuri .
  ?obs dim:populationProjectionVariant ?projectionuri .
	?obs measures:count ?nPeople .
  
	?areauri rdfs:label ?area .
	?perioduri rdfs:label ?periodname .
  	?ageuri rdfs:label ?age .
  	?sexuri rdfs:label ?sex .
  	?projectionuri rdfs:label ?projection . 
	FILTER (?periodname = '2020' && 
            regex( str(?areauri), 'S92|S08|S12' ) &&
  			regex( str(?projectionuri), 'principal')).
	}
  "
  
  return(SPARQL(endpoint,query))
}