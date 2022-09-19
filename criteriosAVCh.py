# -*- coding: utf-8 -*-
# pylint: disable=unnecessary-semicolon
from operator import itemgetter
import math
import numpy as np

class Criterio:
    
    def testarCriterio(self, tomografia, regiao):
        return self.executarTeste(tomografia, regiao);
        
    def executarTeste(self, tomografia, regiao):
        pass; 
    
    def getDescricao(self):
        pass;
    
class ResultadoProcessamentoCriterio:
    """ Guarda as informações do processamento de um criterio de processamento de regiao """   
    
    def __init__(self, passou, resultado):
        self.passou = passou;
        self.resultado = resultado;

    def getResultado(self):
        return self.resultado;

class CriterioAVChTamanhoRegiao(Criterio):
    
    TAMANHO_MINIMO = 150; #voxels, minimo de tamanho, corresponde a aproximadamente TAMANHO_MINIMO*0,298 mm³ 
    
    def getDescricao(self):
        Criterio.getDescricao(self);
        return 'Criterio de tamanho de região';
    
    def executarTeste(self, tomografia, regiao):
        Criterio.executarTeste(self, tomografia, regiao);
        passou = True;
        resultado = "Região atende ao tamanho mínimo necessário(%d) = %d" % (self.TAMANHO_MINIMO, len(regiao.pontos));
        # Análise de tamanho de cada regiao
        if len(regiao.pontos) < self.TAMANHO_MINIMO:
            passou = False;
            resultado = "Regiao não atende ao tamanho mínimo necessário(%d) = %d" % (self.TAMANHO_MINIMO, len(regiao.pontos));
            
        return ResultadoProcessamentoCriterio(passou, resultado);
    
class CriterioAVChDensidadeMinima(Criterio):
    
    DENSIDADE_MINIMA = 46.85;
    
    def getDescricao(self):
        Criterio.getDescricao(self);
        return 'Criterio de densidade mínima';
    
    def executarTeste(self, tomografia, regiao):
        Criterio.executarTeste(self, tomografia, regiao);
        passou = True;
        resultado = "Região atende a densidade mínima necessária %.2f > (densidade mínima = %.2f)" % (regiao.densidadeRegiao, self.DENSIDADE_MINIMA);
        # Análise da densidade de cada região
        if regiao.densidadeRegiao < self.DENSIDADE_MINIMA:
            passou = False;
            resultado = "Região não atende a densidade mínima necessária %.2f < (densidade mínima = %.2f)" % (regiao.densidadeRegiao, self.DENSIDADE_MINIMA);
            
        return ResultadoProcessamentoCriterio(passou, resultado);    
    
class CriterioAVChLimiteDimensional(Criterio):
    
    LIMITE_MINIMO_Z = 3;
    LIMITE_MAXIMO_Z = 110; 
    LIMITE_MINIMO_X = 7;
    
    def getDescricao(self):
        Criterio.getDescricao(self);
        return 'Criterio de limites dimensionais de uma região';
    
    def executarTeste(self, tomografia, regiao):
        Criterio.executarTeste(self, tomografia, regiao);
        passou = True;
        resultado = "Região atende ao critério de dimensão Z %d <= %d <= %d<br/>" % (self.LIMITE_MINIMO_Z, regiao.dimensaoZ, self.LIMITE_MAXIMO_Z);
        resultado += "Região atende ao critério de dimensão X %d <= %d" % (self.LIMITE_MINIMO_X, regiao.dimensaoX);
        # Análise de tamanho de cada regiao
        if regiao.dimensaoZ < self.LIMITE_MINIMO_Z or regiao.dimensaoZ > self.LIMITE_MAXIMO_Z or regiao.dimensaoX < self.LIMITE_MINIMO_X:
            passou = False;
            resultado = "Região não atende ao critério de dimensão Z %d <= %d <= %d<br/>" % (self.LIMITE_MINIMO_Z, regiao.dimensaoZ, self.LIMITE_MAXIMO_Z);
            resultado += "Região não atende ao critério de dimensão X %d < %d" % (regiao.dimensaoX, self.LIMITE_MINIMO_X);
            
        return ResultadoProcessamentoCriterio(passou, resultado);    
    
class CriterioAVChPertencimentoCaixaCraniana(Criterio):

    QTD_MINIMA_PONTOS_CRANIANOS = 3;
    
    def getDescricao(self):
        Criterio.getDescricao(self);
        return 'Criterio de pertencimento a caixa craniana';       
    
    def executarTeste(self, tomografia, regiao):
        self.pontosAnalisados = "";
        self.pontosOsso = 0;        
        Criterio.executarTeste(self, tomografia, regiao);
        passou = True;
        xMedioRegiao = regiao.centroRegiaoX;
        yMedioRegiao = regiao.centroRegiaoY;
        zMedioRegiao = regiao.centroRegiaoZ;
        xMedioTomografia = 255;
        yMedioTomografia = 255;
        zMedioTomografia = int((tomografia.getUltimoCorte()-tomografia.getPrimeiroCorte())/2);

        #Procuro por ao menos 3 pontos com densidade de osso ou AR no caminho 
        # do centro da região ao centro da imagem
        for passo in range(0, 101,2): #de 0 a 100 com passos de 2 (50 iteracoes)
            x = int(xMedioRegiao+(passo/100.0*(xMedioTomografia-xMedioRegiao)));
            y = int(yMedioRegiao+(passo/100.0*(yMedioTomografia-yMedioRegiao)));
            z = int(zMedioRegiao+(passo/100.0*(zMedioTomografia-zMedioRegiao)));
            densidade = tomografia.cortesOriginais[z][y][x]-1024;
            self.pontosAnalisados += "(%d,%d,%d) = %d, " % (x,y,z+tomografia.getPrimeiroCorte()+1,densidade);
            if densidade >= 150 or densidade <= -15:
                self.pontosOsso += 1;

            if self.pontosOsso == self.QTD_MINIMA_PONTOS_CRANIANOS:
                passou = False;
                break;
        
        if passou:
            return ResultadoProcessamentoCriterio(passou, "Região pertence a caixa craniana. "+self.pontosAnalisados);
        else:
            return ResultadoProcessamentoCriterio(passou, "Região não pertence a caixa craniana. "+self.pontosAnalisados);
       
class CriterioAVChRelacaoAreaContentora(Criterio):
    
    MEDIA_CORTES_MINIMA = 0.15;
    
    def getDescricao(self):
        Criterio.getDescricao(self);
        return 'Criterio de relação entre as áreas que contém a regiao e a regiao, por corte';
    
    def executarTeste(self, tomografia, regiao):
        Criterio.executarTeste(self, tomografia, regiao);
        
        #Em um loop por todos os cortes, eu calculo a relação entre a area ocupada pela regiao
        #e o menor retangulo que a contém. Faço a média de cada corte e depois uma média geral
        #dos cortes. Uma área que contém sangue tende a ter um valor alto nesse critério
        mediasCortes = 0;
        listaMedias = [];
        for z in range(regiao.primeiroPontoZ[0], regiao.ultimoPontoZ[0]+1):
            pontosCorte = list(filter (lambda l : l[0] == z , regiao.pontos));
            #Ordenando os pixels de cada corte
            regiaoOrdenadaY = sorted(pontosCorte, key=itemgetter(1));
            regiaoOrdenadaX = sorted(pontosCorte, key=itemgetter(2));
    
            #Adquirindo as dimensões da região(o menor retangulo que contem toda a area desse corte)
            primeiroPontoY, ultimoPontoY = regiaoOrdenadaY[0], regiaoOrdenadaY[len(regiaoOrdenadaY)-1];
            dimensaoY = ultimoPontoY[1]-primeiroPontoY[1]+1;        
            primeiroPontoX, ultimoPontoX = regiaoOrdenadaX[0], regiaoOrdenadaX[len(regiaoOrdenadaX)-1];
            dimensaoX = ultimoPontoX[2]-primeiroPontoX[2]+1;
              
            relacaoArea = 1.0*len(pontosCorte)/(dimensaoX*dimensaoY);
            if relacaoArea < self.MEDIA_CORTES_MINIMA:
                #se a relação for menor, eu apago esse corte da regiao
                regiao.pontos = list(filter (lambda l : l[0] != z , regiao.pontos));
                continue;
            mediasCortes += relacaoArea/regiao.dimensaoZ;
            listaMedias.append("%.2f" % (relacaoArea));
        
        if mediasCortes >= self.MEDIA_CORTES_MINIMA:
            return ResultadoProcessamentoCriterio(True, "Média da relação conteudo/contenedor adequada %.2f >= %.2fF<br/>    *** lista das medias de relação de conteúdo %s" % (mediasCortes, self.MEDIA_CORTES_MINIMA, listaMedias[1:-1]));
        
        else:
            return ResultadoProcessamentoCriterio(False, "Média da relação conteudo/contenedor inadequada %.2f < %.2fF" % (mediasCortes, self.MEDIA_CORTES_MINIMA));

class CriterioAVChEntornoNegativo(Criterio):
    
    BORDA_ENTORNO = 5;
    LIMITE_ENTORNO_NEGATIVO = -15;
    LIMIAR_CONTRASTE = 0.83;
    LIMIAR_RELACAO_PONTOS_OSSO = 1;
    LIMIAR_RAZAO_PONTOS_NEGATIVOS = 0.006;
    mensagem = "";
    
    def getDescricao(self):
        Criterio.getDescricao(self);
        return 'Criterio de entorno negativo e contraste de entorno';    
    
    def testarEntornoNegativo(self, tomografia, regiao, x, y, z):
        for xEntorno in range (x-self.BORDA_ENTORNO, x+self.BORDA_ENTORNO+1):
            for yEntorno in range(y-self.BORDA_ENTORNO, y+self.BORDA_ENTORNO+1):
                for zEntorno in range(z-2, z+3):
                    regiaoDoEntorno = tomografia.regiaoPorPixel.get((zEntorno, yEntorno, xEntorno), -1);
                    if ( (zEntorno, yEntorno, xEntorno) not in regiao.pontosEntorno and 
                        regiaoDoEntorno != regiao.numero and (
                        zEntorno >= 0 and yEntorno >= 0 and xEntorno >= 0) and (
                            zEntorno < len(regiao.cortesOriginais) )) :
                        try:
                            densidadePontoEntorno = regiao.cortesOriginais[zEntorno][yEntorno][xEntorno]-1024;
                            #print("Testando o ponto (%d, %d, %d) = %d" % (zEntorno+regiao.primeiroCorte+1, yEntorno, xEntorno, densidadePontoEntorno));
                            #Adiciono o ponton que está sendo testado a lista de pontos de entorno para esta regiao
                            regiao.pontosEntorno[(zEntorno, yEntorno, xEntorno)] = densidadePontoEntorno;
                            if not densidadePontoEntorno > 100 and not densidadePontoEntorno < self.LIMITE_ENTORNO_NEGATIVO:
                                self.densidadesEntorno.append(densidadePontoEntorno);
                            self.listaDensidadesPontosEntorno.append("(%d,%d,%d) = %d" % ( xEntorno, yEntorno, zEntorno+regiao.primeiroCorte+1, densidadePontoEntorno));
                            if densidadePontoEntorno > 100:
                                self.pontosDeOsso.append((zEntorno, yEntorno, xEntorno));
                            #Se houver algum ponto com densidade negativa, signifca que a area de entorno
                            #contém algum ponto com Ar, o que não deve acontecer com uma regiao intraparenquimatosa
                            elif densidadePontoEntorno < self.LIMITE_ENTORNO_NEGATIVO and z == zEntorno:
                                self.densidadesPontosNegativos.append(densidadePontoEntorno);
                                #self.mensagem += "Encontrado um ponto(%d, %d, %d) do entorno de (%d, %d, %d) negativo: %d<br/>" % (xEntorno, yEntorno, zEntorno+regiao.primeiroCorte+1, x, y, z+regiao.primeiroCorte+1, densidadePontoEntorno);
                            else:
                                self.densidadeRegiaoEntorno += densidadePontoEntorno;
                        except IndexError:
                            densidadePontoEntorno = 0;
                    
        return False;    
    
    def executarTeste(self, tomografia, regiao):
        Criterio.executarTeste(self, tomografia, regiao);
        self.entornoNegativo = False;
        self.densidadeRegiaoEntorno = 0;   
        self.densidadesPontosNegativos = [];   
        self.pontosDeOsso = [];  
        self.densidadesEntorno = [];
        self.mensagem = "";
        self.razaoPontosNegativos = 0;
                
        #calculos a partir dos valores de cada voxel
        self.listaDensidadesPontosEntorno = [];
        for ponto in regiao.pontos:
            if not self.entornoNegativo:
                self.entornoNegativo = self.testarEntornoNegativo(tomografia, regiao, ponto[2], ponto[1], ponto[0]);
        self.volumeRegiaoEntorno = len(regiao.pontosEntorno); 
        self.densidadeRegiaoEntorno /= 1.0*self.volumeRegiaoEntorno+1;
        self.relacaoPontosOsso = 1.0*len(regiao.pontos)/(len(self.pontosDeOsso)*1.0+0.00000001);

        dpe = np.array(self.densidadesEntorno);
        self.medianaEntorno = np.median(dpe);
        self.desvioPadraoEntorno = np.std(dpe);     

        self.mensagem += "Descrição do Entorno:<br/>";
        self.mensagem += "Densidade do entorno %.2f<br/>" % (self.densidadeRegiaoEntorno);   
        self.mensagem += "Mediana %.2f - Desvio Padrão %.5f<br/>" % (self.medianaEntorno, self.desvioPadraoEntorno);
        if len(self.densidadesPontosNegativos) == 0:
            self.mensagem += "Não há pontos de densidade negativa no entorno da região.<br/>";
        else:
            self.razaoPontosNegativos = 1.0*len(self.densidadesPontosNegativos)/self.volumeRegiaoEntorno;
            self.mensagem +=  "Existem %d pontos de densidade negativa no entorno da região. Quantidade total de pontos / pontos de densidade negativa = %f.<br/>" % (len(self.densidadesPontosNegativos), self.razaoPontosNegativos);
            if self.razaoPontosNegativos > self.LIMIAR_RAZAO_PONTOS_NEGATIVOS:
                self.entornoNegativo = True;
                
        #Testa a relação de contraste entre a regiao o entorno
        #deve haver um contraste suficiente entre a regiao e o seu entorno
        #se o contraste for suficiente, entao significa que visivelmente, um humano
        #consegue distinguir uma região do seu entorno
        self.contraste = False;
        self.ultrapassouLimiarOsso = False;
        if self.densidadeRegiaoEntorno <= self.LIMIAR_CONTRASTE*regiao.densidadeRegiao:
            self.mensagem += "Contraste suficiente entre a região e seu entorno %.2f (%.2f%%) <= %.2f%%<br/>" % (self.densidadeRegiaoEntorno, self.densidadeRegiaoEntorno/regiao.densidadeRegiao, self.LIMIAR_CONTRASTE);            
            #self.mensagem += "    Lista de densidades(%d) : %s" % (len(self.listaDensidadesPontosEntorno), ', '.join(self.listaDensidadesPontosEntorno));                        
            self.contraste = True;
        else:
            self.mensagem += "Contraste da região de entorno insuficiente %.2f (%.2f%%) > %.2f%%<br/>" % (self.densidadeRegiaoEntorno, self.densidadeRegiaoEntorno/regiao.densidadeRegiao*100, self.LIMIAR_CONTRASTE);
            #self.mensagem += "    Lista de densidades(%d) : %s" % (len(self.listaDensidadesPontosEntorno), ', '.join(self.listaDensidadesPontosEntorno));            

        #Uma regiao intraparenquimatosa não pode conter um ponto de entorno muito negativo
        if not self.entornoNegativo and self.volumeRegiaoEntorno != 0:
            self.mensagem += "Razão de pontos de densidade negativa menor que o limiar (%f < %f). Densidade regiao entorno = %.2f<br/><br/>" % (self.razaoPontosNegativos, self.LIMIAR_RAZAO_PONTOS_NEGATIVOS, self.densidadeRegiaoEntorno);
        else:
            self.mensagem += "Razão de pontos de densidade negativa MAIOR que o limiar (%f > %f). Densidade regiao entorno = %.2f<br/><br/>" % (self.razaoPontosNegativos, self.LIMIAR_RAZAO_PONTOS_NEGATIVOS, self.densidadeRegiaoEntorno);

        #A relação entre os pontos da regiao e o de osso deve ser que o limiar estabelecido
        #essa relação estabelece que a quantidade de pontos da região é muito maior do que a quantidade
        #de pontos com densidade semelhante a osso no entorno da região, o que caracteriza uma
        #região intraparenquimatosa
        if self.relacaoPontosOsso >= self.LIMIAR_RELACAO_PONTOS_OSSO:
            self.mensagem += "Relação pontos osso maior que o limiar (%f > %f)<br/>" % (self.relacaoPontosOsso, self.LIMIAR_RELACAO_PONTOS_OSSO);            
        else:
            self.mensagem += "Relação pontos osso menor que o limiar (%f < %f)<br/>" % (self.relacaoPontosOsso, self.LIMIAR_RELACAO_PONTOS_OSSO);                        
            self.ultrapassouLimiarOsso = True;

        return ResultadoProcessamentoCriterio(not self.entornoNegativo and self.contraste and not self.ultrapassouLimiarOsso,  self.mensagem);   
        
class CriterioCircularidadePequenasAreas(Criterio):
    
    TAMANHO_MAXIMO = 100;

    def caminharVetor(self, pInicial, angulacaoVetor):
        proximoPonto = pInicial;
        (y,x) = (proximoPonto[0], proximoPonto[1]);
        if abs(angulacaoVetor) % 360 != 0 and abs(angulacaoVetor) % 360 !=  180:
            proximoPonto = (y+1 if (angulacaoVetor+360) % 360 < 180 and (angulacaoVetor+360) % 360 > 0 else y-1, x);
        if abs(angulacaoVetor) % 360 != 90 and abs(angulacaoVetor) % 360 != 270:
            proximoPonto = (proximoPonto[0], x+1 if (angulacaoVetor+360) % 360 < 90 or (angulacaoVetor+360) % 360 > 270 else x-1);
        return proximoPonto;    
    
    def getDescricao(self):
        Criterio.getDescricao(self);
        return 'Criterio de circularidade de pequenas areas';
    
    def executarTeste(self, tomografia, regiao):
        Criterio.executarTeste(self, tomografia, regiao);
        #O proposito desse teste é verificar se, caso seja uma area pequena, um teste de circularidade pode
        #revelar que se trata de uma regiao que não corresponde ao espalhamento habitual de um liquido
        
        listaCircularidades = [];
        for z in range(regiao.primeiroPontoZ[0], regiao.ultimoPontoZ[0]+1):
            pontosCorte = list(filter (lambda l : l[0] == z , regiao.pontos));
            area = len(pontosCorte);
            if area > self.TAMANHO_MAXIMO or area == 0:
                continue;
            
            perimetro = 1;
            angulacaoVetor = 0; #angulacao em graus
            #Nesse Momento eu tento percorrer o contorno da regiao, para calcular seu perimetro
            #O objetivo é percorrer o contorno no sentido horario
            #começo do ponto inicial, e termino quando chegar nele novamente
            
            ultimoPonto = pInicial = pontosCorte[0];
            regiaoContorno = [];
            regiaoContorno.append(pInicial);
            proximoPonto = self.caminharVetor(ultimoPonto[1:3], angulacaoVetor);
            while ultimoPonto != pInicial or perimetro == 1:
                if tomografia.regiaoPorPixel.get((z, proximoPonto[0], proximoPonto[1]), -1) == regiao.numero:
                    regiaoContorno.append(proximoPonto);
                    perimetro += 1;
                    if perimetro >= area:
                        break;
                    angulacaoVetor -= 90;
                    ultimoPonto = (z, proximoPonto[0], proximoPonto[1]);
                else:
                    angulacaoVetor += 45;
                    if angulacaoVetor >= 360:
                        print("*** Erro no calculo do perimetro da regiao %d, corte %d" % (regiao.numero, z));                        
                        break;
                proximoPonto = self.caminharVetor(ultimoPonto[1:3], angulacaoVetor);
                
            #Fórmula da circularidade
            circularidade = 4.0 * math.pi * area / perimetro**2;
            if circularidade > 1:
                regiao.pontos = [p for p in regiao.pontos if p not in pontosCorte];
            
            listaCircularidades.append(circularidade);
        
        #if mediasCortes >= self.MEDIA_CORTES_MINIMA:
            #return ResultadoProcessamentoCriterio(True, "Média da relação conteudo/contenedor adequada %f >= %fF<br/>    *** lista das medias %s" % (mediasCortes, self.MEDIA_CORTES_MINIMA, str(listaMedias)[1:-1]));
        #else:
        return ResultadoProcessamentoCriterio(True, "Circularidades *** %s" % (str(listaCircularidades)[1:-1]));
        


 