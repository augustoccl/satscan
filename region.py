def Enum(*sequential, **named):
    Enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), Enums)
    
TipoRegiao = Enum('TipoRegiao', 'INDETERMINADA', 'SANGUE', 'PSEUDO_SANGUE');

class Region:
        
    def __init__(self, numero, tipo, qtdCortes, cortesOriginais, primeiroCorte):
        self.tipo = tipo;
        self.numero = numero;
        self.cortesOriginais = cortesOriginais;
        self.pontos = [];
        self.pontosEntorno = {};
        self.qtdCortes = qtdCortes;
        self.primeiroCorte = primeiroCorte;
        self.resultadosProcessamentos = [];
        
    def getDescricao(self):
        descricao  = "*** Regiao %d - Volume %.0f - Densidade %.2f<br/>" % (self.numero, self.volumeRegiao, self.densidadeRegiao);
        descricao += "    Mediana = %d || Moda = %d<br/>" % (self.medianaRegiao, self.modaRegiao);
        descricao += "    Desvio Padrão = %.5f<br/>" % (self.desvioPadrao);
        descricao += "    DimensãoX = %d DimensãoY = %d DimensãoZ = %d<br/>" % (self.dimensaoX, self.dimensaoY, self.dimensaoZ);
        descricao += "    PrimeiroX = %d UltimoX = %d<br/>" % (self.primeiroPontoX[2], self.ultimoPontoX[2]);
        descricao += "    PrimeiroY = %d UltimoY = %d<br/>" % (self.primeiroPontoY[1], self.ultimoPontoY[1]);
        descricao += "    PrimeiroZ = %d UltimoZ = %d<br/>" % (self.primeiroPontoZ[0]+self.primeiroCorte+1, self.ultimoPontoZ[0]+self.primeiroCorte+1);
        descricao += "    MenorX = %d MaiorX = %d<br/>" % (self.menorX, self.maiorX);
        descricao += "    InicialY = %d FinalY = %d<br/>" % (self.inicialY, self.finalY);
        descricao += "    InicialZ = %d FinalZ = %d<br/>" % (self.inicialZ, self.finalZ);
        descricao += "    --------------------------<br/>";
        
        return descricao; 

    def adicionarResultado(self, rp):
        self.resultadosProcessamentos.append(rp);

    def getDescricaoComResultados(self):
        r = self.getDescricao();
        for rp in self.resultadosProcessamentos:
            if rp.passou:
                r += rp.resultado+'<br/>';
            else:
                r += f"""<span style='color:red'>{rp.resultado}</span><br/>""";
        return r;

    def prepararProcessamento(self):
        #Ordenando os pixels de cada regiao
        regiaoOrdenadaZ = sorted(self.pontos, key=itemgetter(0));
        regiaoOrdenadaY = sorted(self.pontos, key=itemgetter(1));
        regiaoOrdenadaX = sorted(self.pontos, key=itemgetter(2));

        #Adquirindo as dimensões da região(o menor cubo que contem toda a região)
        self.primeiroPontoZ, self.ultimoPontoZ = regiaoOrdenadaZ[0], regiaoOrdenadaZ[len(regiaoOrdenadaZ)-1];
        self.dimensaoZ = self.ultimoPontoZ[0]-self.primeiroPontoZ[0]+1;        
        self.primeiroPontoY, self.ultimoPontoY = regiaoOrdenadaY[0], regiaoOrdenadaY[len(regiaoOrdenadaY)-1];
        self.dimensaoY = self.ultimoPontoY[1]-self.primeiroPontoY[1]+1;        
        self.primeiroPontoX, self.ultimoPontoX = regiaoOrdenadaX[0], regiaoOrdenadaX[len(regiaoOrdenadaX)-1];
        self.dimensaoX = self.ultimoPontoX[2]-self.primeiroPontoX[2]+1;  

        self.centroRegiaoZ = int(self.primeiroPontoZ[0]+self.dimensaoZ/2);
        self.centroRegiaoY = int(self.primeiroPontoY[1]+self.dimensaoY/2);
        self.centroRegiaoX = int(self.primeiroPontoX[2]+self.dimensaoX/2);    

        self.centroImagemX, self.centroImagemY, self.centroImagemZ = (255,255, self.qtdCortes/2);
        #Primeiro eu encontro o quadrante que o centro dessa regiao se encontra em relação ao centro da imagem
        #Valores iniciais considerando o quadrante superior esquerdo 
        self.menorX, self.maiorX = (self.centroRegiaoX, self.centroImagemX);
        self.inicialY, self.finalY = (self.centroRegiaoY, self.centroImagemY);
        self.inicialZ, self.finalZ = self.centroRegiaoZ, self.centroImagemZ;
        if self.centroRegiaoX > self.centroImagemX:
            self.menorX, self.maiorX = self.centroImagemX, self.centroRegiaoX;       
            self.inicialY, self.finalY = self.centroImagemY, self.centroRegiaoY;
            self.finalZ, self.inicialZ = self.centroImagemZ, self.centroRegiaoZ;

        self.densidadeRegiao = 0;
        self.volumeRegiao = 0;
        self.densidadePontos = [];
        for ponto in self.pontos:
            self.volumeRegiao += 1;
            d = self.cortesOriginais[ponto[0]][ponto[1]][ponto[2]]-1024;
            self.densidadeRegiao += d;
            self.densidadePontos.append(d);
        dp = np.array(self.densidadePontos);
        self.densidadeRegiao /= 1.0*self.volumeRegiao;
        self.medianaRegiao = np.median(dp);
        self.modaRegiao = np.apply_along_axis(lambda x: np.bincount(x).argmax(), axis=0, arr=dp);
        self.desvioPadrao = np.std(dp);