



## Las redes transcripcionales son redes dirigidas donde los nodos representan genes 
## y se traza una arista de un nodo g1 al node g2 cuando g1 codifica por un factor de 
## transcripci?n se une al promoter del g2. T?picamente, est?s redes se construyen a 
## partir de datos de ChIP-seq y suponen un refinameinto de la redes de co-expresi?n 
## g?nica. 

## Este script recibe como entrada un fichero donde a su vez se guardan los 
## identificadores de los genes para los cuales a partir de datos de ChIP-seq se han 
## determinado sus dianas potenciales.

tf.name <- read.table(file = "target.txt", header = FALSE)[[1]]
length(tf.name)

## Para obtener el conjunto total de genes de la red leemos todos los ficheros 
## correspondientes a la dianas potenciales de los factores de transcripci?n y los
## vamos almacenando en un vector. 

genes.in.network <- c()
for(i in 1:length(tf.name))
{
  current.tf.file <- paste(tf.name[i],".txt",sep="")
  current.tf.data <- read.table(file=current.tf.file,as.is=TRUE)
  current.tf.targets <- current.tf.data[,1]
  genes.in.network <- c(genes.in.network,current.tf.targets)
}

genes.in.network <- unique(genes.in.network)
genes.in.network
length(genes.in.network)

## Generamos incialmente una matriz de adyacencia con todos los valores 0 y vamos
## rellenándola: los valores 1 son las dianas de los correspondientes factores
## de transcripci?n. 
ad.matrix<- matrix(data=0, 
                   nrow = length(genes.in.network), 
                   ncol = length(genes.in.network))

rownames(ad.matrix)<- genes.in.network
colnames(ad.matrix) <- genes.in.network

ad.matrix[1:7,1:7]

for(i in 1:length(tf.name))
{
  current.tf <- tf.name[i]
  current.tf.file <- paste(current.tf, ".txt",sep="")
  current.tf.data <- read.table(file=current.tf.file,as.is=TRUE)
  current.tf.targets <- current.tf.data[,1]
  ad.matrix[current.tf,current.tf.targets] <- 1
}

## Finalmente, construimos la red transcripcional a partir de su matriz de adyacencia
## y la guardamos en formato gml. 

library(igraph)

gene.transcriptional.network <- graph.adjacency(ad.matrix, mode="directed")
write.graph(gene.transcriptional.network,file="transcriptional_gene_network.gml",format="gml")

gene.transcriptional.network

##AQUÍ SE GENERA LA RED TRANSCRIPCIONAL PERO NO ES LA QUE HAY QUE USAR EN EL CYTOSCAPE

## Esta red es excesivamente grande y por motivos de recursos computacionales y tiempo
## vamos a trabajar s?lo con la red inducida que contiene los factores de transcripci?n
## y las aristas entre ellos.

##AQUÍ SÍ SE GENERA LA RED TRANSCRIPCIONAL QUE HAY QUE USAR EN EL CYTOSCAPE.
##CUANDO ESO ESTÉ SI QUIERES ME LO PASAS Y LO HAGO A VER SI ME VA

tfs.network <- induced.subgraph(gene.transcriptional.network,tf.name)
plot.igraph(x = tfs.network,vertex.size=8,edge.arrow.size=0.5,vertex.label="",vertex.color="blue")
write.graph(tfs.network,file="red_transcripcional_tf_inducida.gml",format="gml")

tfs.network<- read.graph(file = "red_transcripcional_tf_inducida.gml", format = "gml" )

## Uno de los primeros an?lisis que se realizan sobre redes transcripcionales es la 
## búqueda de subgrafos o patrones no aleatorios. Estos subgrafos se denominan 
## motivos de red.
## Formalmente, un motivo de red es un subgrafo que aparece un n?mero significativamente
## de veces mayor en la red de interés que en redes aleatorias que cumplen las 
## mismas propiedades.

## Primero, realizamos un an?lisis r?pido para comprobar que esta red inducida
## NO es libre de escala. 

# Calculo del grado de los nodos

network.degrees <- degree(tfs.network)
hist(network.degrees,col="blue",xlab="Node degree", ylab="Probability",main="Degree distribution")

# Calculo de la frecuencia absoluta del grado de los nodos
degree.frequencies <- table(network.degrees)
# Eliminamos nodos con grado 0 para poder aplicar log10

# Transformaci?n logar?tmica
log10.degrees.frequencies <- log10(degree.frequencies)
log10.node.degrees <- log10(as.numeric(names(degree.frequencies)))

# Regresi?n lineal
lm.r <- lm(log10.degrees.frequencies ~ log10.node.degrees)
summary(lm.r)

##EN EL SUMMARY HAY QUE VER EL VALOR DEL P-VALOR. Si es menor que 0.05 es libre de escala
## El p-valor sale 0.9898, NO es libre de escala

## Las redes aleatorias que cumplen las mismas propiedas que nuestra
## red de interés son las generadas seg?n el modelo de Erdos-Renyi. La funci?n
## erdos.renyi.game genera redes aleatorias que siguen una distribuci?n de Poisson
## con un n?mero de nodos y aristas dado. 
tfs.network

##AQUÍ TIENES QUE METER EL NÚMERO DE NODOS EN EL n (EL PRIMER NÚMERO QUE SALE)

random.graph <- erdos.renyi.game(n=10, p.or.m=22, type="gnm", directed=TRUE, loops=TRUE)
plot.igraph(x = random.graph,vertex.size=8,edge.arrow.size=0.5,vertex.label="",vertex.color="blue")

## Empezamos viendo si es un motivo de red la autorregulaci?n.
## Para ver el n?mero de genes autorregulados en una red basta con sumar los
## elementos de la diagonal principal. 

## N?mero de genes autorregulados en la red de factores de transcripci?n. 
tfs.adjacency <- as.matrix(get.adjacency(tfs.network))
tfs.adjacency
diag(tfs.adjacency)
autorregulation.in.tfs <- sum(diag(tfs.adjacency))
autorregulation.in.tfs

## N?mero de genes autorregulados en la red 
autorregulation.in.random <- sum(diag(as.matrix(get.adjacency(random.graph))))
autorregulation.in.random ##•ESTE ES UN NÚMERO QUE SE GENERA DE FORMA ALEATORIA

## Para comprobar la signficancia de la autorregulaci?n generamos 1000 redes
## aleatorias y calculamos el n?mero de autorregulaciones en cada una almacen?ndola
## en un vector.
autorregulation.random.graphs <- vector(length=10000, mode="numeric")

for (i in 1:10000)
{
  random.graph <- erdos.renyi.game(n=21, p.or.m=115, type="gnm", directed=TRUE, loops=TRUE)
  autorregulation.random.graphs[i] <- sum(diag(as.matrix(get.adjacency(random.graph))))
}

mean(autorregulation.random.graphs)
sd(autorregulation.random.graphs)

sum(autorregulation.random.graphs > autorregulation.in.tfs)/10000 ##El p.valor es de 10^-4 que es el numero de veces que he generado

##AQUÍ AL FINAL LO QUE SE MIRA ES LA RELEVANCIA DE LOS MOTIVOS AUTORREGULADOS QUE DEPENDE
##DEL VALOR QUE SALGA EN LAS REDES GENERADAS AUTOMÁTICAMENTE


## Motivos de red de tres nodos

## La funci?n graph.motifs recibe como entrada una red y un tama?o de subgrafo k
## (en la actualidad s?lo puede recibir tama?os 3 o 4) y devuelve el n?mero de
## veces que se encuentra cada subgrafo con k nodos en la red. 

occurrency.subgraph.three.nodes <- graph.motifs(tfs.network, size=3)
occurrency.subgraph.three.nodes
## Determina cuantos tipos de subgrafos hay; NA es que no aparece ninguna vez

length(occurrency.subgraph.three.nodes)

## La funci?n graph.isocreate genera todos los grafos posibles de un tama?o dado
## empezando a enumerarlos por 0. 

##LOS GRAFOS SON LOS DIFERENTES MOTIVOS QUE DIMOS EL AÑO PASADO (EN PLAN QUE UN REGULA A OTRO Y TAL)
##ENTONCES EXISTEN 16 DIFERENTES Y NOSOTROS TENEMOS QUE MIRAR LOS QUE SALEN (EN PLAN SI SALE NA ES QUE NO 
##APARECE NINGUNA VEZ)

plot.igraph(graph.isocreate(size=3, number=4))
plot.igraph(graph.isocreate(size=3, number=6))
plot.igraph(graph.isocreate(size=3, number=7))
plot.igraph(graph.isocreate(size=3, number=8))
plot.igraph(graph.isocreate(size=3, number=9))
plot.igraph(graph.isocreate(size=3, number=12))
plot.igraph(graph.isocreate(size=3, number=13))



## Para ver si cada uno de los 16 posibles subgrafos es o no motivo de red
## realizamos le procedimiento anterior basado en la generaci?n de redes aleatorias.
## En este caso el acumulador es matricial. 
motifs.3.random.graph <- matrix(0,nrow=1000, ncol=16)
motifs.3.random.graph[1:3,]

##AQUÍ CAMBIA EL N
for (i in 1:1000)
{
  random.graph <- erdos.renyi.game(n=10, p.or.m=22, type="gnm", directed=TRUE, loops=TRUE)
  motifs.3.random.graph[i,] <- graph.motifs(random.graph, size=3)
}

motifs.3.random.graph[1:3,]

##ESTAMOS GENERANDO REDES ALEATIRUAS PARAVER SI ES MOTIVO DE RED

## subgrafo 5
plot(graph.isocreate(size=3, number=4))
occurrency.subgraph.three.nodes[5]

mean(motifs.3.random.graph[,5])
sd(motifs.3.random.graph[,5])
sum(motifs.3.random.graph[,5] >1)/1000

##No significativo pq el p.valor que yo he generado es de 10^-4 (en el ejemplo de fran)

## subgrafo 7
plot(graph.isocreate(size=3, number=6))
occurrency.subgraph.three.nodes[7]

mean(motifs.3.random.graph[,7])
sd(motifs.3.random.graph[,7])
sum(motifs.3.random.graph[,7] > 42)/1000 ##0

## Subgrafo 8 
plot(graph.isocreate(size=3, number=7))
occurrency.subgraph.three.nodes[8]

mean(motifs.3.random.graph[,8])
sd(motifs.3.random.graph[,8])
sum(motifs.3.random.graph[,8] > 8)/1000

## Subgrafo 9
plot(graph.isocreate(size=3, number=8))
occurrency.subgraph.three.nodes[9]

mean(motifs.3.random.graph[,9])
sd(motifs.3.random.graph[,9])
sum(motifs.3.random.graph[,9] > 1)/1000

## Subgrafo 10
plot(graph.isocreate(size=3, number=9))
occurrency.subgraph.three.nodes[10]

mean(motifs.3.random.graph[,10])
sd(motifs.3.random.graph[,10])
sum(motifs.3.random.graph[,10] > 14)/1000 ##0

## Subgrafo 13 
plot(graph.isocreate(size=3, number=12))
occurrency.subgraph.three.nodes[13]

mean(motifs.3.random.graph[,13])
sd(motifs.3.random.graph[,13])
sum(motifs.3.random.graph[,13] > 1)/1000

## Subgrafo 14
plot(graph.isocreate(size=3, number=13))
occurrency.subgraph.three.nodes[14]

mean(motifs.3.random.graph[,14])
sd(motifs.3.random.graph[,14])
sum(motifs.3.random.graph[,14] > 1)/1000

## Los subgrafos que son motivos con interés biológico son el 7 y el 10

subgrafo.7 <- graph.isocreate(size=3, number=6)

## Para encontrar todas las instancias de este subgrafo en nuestra red usamos 
## la siguiente instrucci?n. 
maps.subgrafo.7 <- graph.subisomorphic.lad(subgrafo.7, tfs.network, all.maps=TRUE)[["maps"]]
maps.subgrafo.7

##ESTO DICE FERNANDO QUE ES QUE TE DICE DÓNDE ESTÁN LOS SUBGRAFOS EN LA RED

subgrafo.10<- graph.isocreate(size = 3, number = 9)
maps.subgrafo.10<- graph.subisomorphic.lad(subgrafo.10, tfs.network, all.maps = TRUE)[["maps"]]
maps.subgrafo.10
