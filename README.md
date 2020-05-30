# rENA
Epistemic Network Analysis (ENA)
For Python developers

# Usando o rENA no python
O exemplo abaixo mostra como utilizar a biblioteca do R rENA com a liguagem python. Todas as bibliotecas do python podem ser usadas para mnipular e arrumar os dados para serem plotados.

Todas as funçoes da biblioteca podem ser chamadas usando a interface rENA, um outro detalhe ,e que os pontos (.) devem ser substituidos por underscore (_).
## Exemplo: 
A função ena.plot() pode ser chamada usando rENA.ena_plot() 
A função ena.plot.points() pode ser chamada usando rENA.ena_plot_points().

# Bibliotecas necessárias:
Devem ser instaladas as bibliotecas R rENA e data.table usando os comandos abaixo:
## install.packages("rENA")
## install.packages("data.table")

No python devemos instalar a biblioteca rpy2 usando o seguinte comando:
pip install rpy2
