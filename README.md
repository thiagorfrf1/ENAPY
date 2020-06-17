# rENA 
<center>
  
![alt text](https://github.com/thiagorfrf1/ENAPY/blob/master/ena.png)

</center>

Epistemic Network Analysis (ENA) <br>

(ENA) é um método para identificar e quantificar conexões entre elementos em dados codificados e representá-los em modelos dinâmicos de rede. Uma característica fundamental da ferramenta ENA é que ela permite que os pesquisadores comparem redes diferentes, visualmente e por meio de estatísticas resumidas que refletem a estrutura ponderada das conexões. A interface também permite que os usuários vejam os dados originais que contribuíram para cada uma das conexões na representação de rede. Assim, a ENA pode ser usada para abordar uma ampla gama de questões de pesquisa qualitativa e quantitativa.

Os pesquisadores usaram a ENA para analisar e visualizar uma ampla gama de fenômenos, incluindo: conexões cognitivas que os alunos fazem enquanto resolvem problemas complexos; interações entre diferentes regiões do cérebro em dados de RMf; coordenação do olhar social; integração de habilidades operatórias durante procedimentos cirúrgicos; e muitos outros.

https://cran.r-project.org/web/packages/rENA/rENA.pdf

# Por que usar a biblioteca?
Uma outra vantagem é que a biblioteca rena original está sendo usada como motor, todas as atualizações e futuras novas funções poderão ser utilizadas com a biblioteca

# Instalando

## Pip:

## pip install enapy

# Usando o rENA no python
O exemplo abaixo mostra como utilizar a biblioteca do R rENA com a liguagem python. Todas as bibliotecas do python podem ser usadas para mnipular e arrumar os dados para serem plotados.

Todas as funçoes da biblioteca podem ser chamadas usando a interface rENA, um outro detalhe ,e que os pontos (.) devem ser substituidos por underscore (_).
## Exemplo: 
A função ena.plot() pode ser chamada usando rENA.ena_plot() <br>
A função ena.plot.points() pode ser chamada usando rENA.ena_plot_points().

# Requisitos:
Python 3.6

rpy2


R 

data.table

rENA

## Bibliotecas necessárias:
Devem ser instaladas as bibliotecas R rENA e data.table usando os comandos abaixo:
## install.packages("rENA")
## install.packages("data.table")

No python devemos instalar a biblioteca rpy2 usando o seguinte comando:
pip install rpy2

## Exemplos
Exemplos estão disponíveis no arquivo ena.ipynb na pasta “examples”.
