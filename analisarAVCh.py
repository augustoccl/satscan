# encoding: utf-8
# pylint: disable=unnecessary-semicolon
'''
satscan.abrirTomografiaDycom -- shortdesc

@author:     Augusto Camarotti

@license:    GPL 3.0

@contact:    augustoccl@gmail.com
'''

import sys
import os
import time
import copy
import png
import numpy as np
import pydicom
import re
import json
from glob import glob
from criteriosAVCh import CriterioAVChTamanhoRegiao;
from criteriosAVCh import CriterioAVChLimiteDimensional;
from criteriosAVCh import CriterioAVChDensidadeMinima;
from criteriosAVCh import CriterioAVChRelacaoAreaContentora;
from criteriosAVCh import CriterioAVChPertencimentoCaixaCraniana;
from criteriosAVCh import CriterioAVChEntornoNegativo;
from resultReports import Reports;
import multiprocessing;
from multiprocessing import Manager, Value, Array;
from multiprocessing.managers import BaseManager;
from functools import partial;
import gc;
import multiprocessing;

from operator import itemgetter, attrgetter
from multiprocessing import freeze_support
from _ast import Or

def Enum(*sequential, **named):
    Enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), Enums)

TipoRegiao = Enum('TipoRegiao', 'INDETERMINADA', 'SANGUE', 'PSEUDO_SANGUE');
diagnosticoTomografia = Enum('SEM_ALTERACOES', 'AVCh');

MINIMA_DENSIDADE_SANGUE = 48;
MAXIMA_DENSIDADE_SANGUE = 78; 
SEQUENCIA_MINIMA = 4;   

class Tomografia:

    cortesOriginais = None;
    regioes = None;
    regiaoPorPixel = None;
    pastaPaciente = None;
    diagnostico = None;
   
    pasta = None;
    
    def getLaudoRadiologista(self):
        return self.laudoRadiologista;

    def setLaudoRadiologista(self, laudo):
        self.laudoRadiologista = laudo;
    
    def getCorteOriginal(self, c):
        return self.cortesOriginais[c];
    
    def getCortesOriginais(self):
        return self.cortesOriginais;
    
    def setPasta(self, p):
        self.pasta = p;
        
    def getNomePaciente(self):
        return self.nomePaciente;
    
    def getPatientID(self):
        return self.patientID;
    
    def getStudyId(self):
        return self.studyId;
    
    def getStudyDate(self):
        return self.studyDate;
    
    def getStudyTime(self):
        return self.studyTime;
    
    def getPrimeiroCorte(self):
        return self.primeiroCorte;

    def getUltimoCorte(self):
        return self.ultimoCorte;

    def getAltura(self):
        return self.altura;
    
    def getLargura(self):
        return self.largura;
    
    def getRegiaoPorPixel(self):
        return self.regiaoPorPixel;
    
    def getRegioes(self):
        return self.regioes;

    def getPastaPaciente(self):
        return self.pastaPaciente;

    def setPaciente(self, pasta):
        self.pastaPaciente = pasta;

    def setRegioes(self, reg):
        self.regioes = reg;
    
    def inicializar(self):
        cortes = [];
        self.cortesOriginais = [];
        self.regiaoPorPixel = {};
        self.totalCortes = 0;
        self.ultimoCorte = 0;
        self.primeiroCorte = 0;
        self.qtdCortes = 0;
        self.cortesOrdenados = None;
        self.cortesModificados = None;
        self.nomePaciente = None;
        self.planilhaPontos = None;
        self.largura = None;
        self.altura = None;
        self.pastaPaciente = "";  
        self.diagnostico = "";      
        #cria uma lista pra guardar as regiões contiguas
        self.regioesContiguas = {};
        
        # Abrir a pasta com as tomografias e ler as imagens
        arquivosDcm = [y for x in os.walk(self.pasta) for y in glob(os.path.join(x[0], '*.dcm'))];
        primeiroDS = None;
        for arquivo in arquivosDcm:
            ds = pydicom.dcmread(arquivo);
            serie = ds.SeriesDescription;
            if serie in ('PARTES MOLES', 'Vol PM', 'PARTE MOLES','PARTE MOLES S/C'):
                cortes.append(ds);
                if (primeiroDS is None):
                    primeiroDS = ds;
                    self.largura, self.altura = primeiroDS.pixel_array.shape;
                    
        if len(cortes) < 1:
            raise ImportError('Nenhum corte de tomografia de partes moles encontrado na pasta');

        cortesOrdenados = sorted(cortes, key=lambda c: c.InstanceNumber);  
        self.totalCortes = len(cortesOrdenados);
        self.ultimoCorte = self.totalCortes*9/10;
        self.primeiroCorte = self.totalCortes/10;
        self.qtdCortes = self.ultimoCorte-self.primeiroCorte;
             
        self.nomePaciente = str(primeiroDS.PatientName);
        self.studyId = str(primeiroDS.StudyID);
        self.studyDate = str(primeiroDS.StudyDate);
        self.studyTime = str(primeiroDS.StudyTime);
        self.patientID = str(primeiroDS.PatientID);
        
        if len(cortesOrdenados) < 1:
            return None;
            
        print ('Total de cortes: %d - Primeiro e Ultimo cortes considerados %d,%d' % (self.totalCortes, self.primeiroCorte, self.ultimoCorte));    
        self.cortesModificados = [];
        for c in range(0, len(cortesOrdenados)):
            if c < self.primeiroCorte or c > self.ultimoCorte:
                continue;
            novaArray = cortesOrdenados[c].pixel_array.tolist();
            self.cortesModificados.append(novaArray);
            self.cortesOriginais.append(copy.deepcopy(novaArray)); 
        print("Tomografia de %s lida com sucesso." % self.nomePaciente);
    
    def selecionarRegioesTomografia(self):
        #===========================================================================
        # Fazendo a selecao de regioes que possivelente contem sangue
        #===========================================================================
        print ("Analisando a tomografia de %s - %s" % (self.nomePaciente, self.laudoRadiologista));
                
        #aqui faremos uma analise em 3 dimensões da tomografia como um todo
        #a variável cortes conterá uma array de cortes
        #suporei que todos os cortes tem a mesma dimensao do primeiro corte
        #a variavel z será utilizada como profundidade, o que equivale a cada corte 
        #IMPORTANTE: Os cortes devem estar em ordem topográfica
        #faco sempre uma copia do corte para nunca trabalhar com o array original
        #volume por pixel
        #cada pixel em um corte tem as dimensões de 0,488 mm² de largura e altura
        #e 1,25mm de profundidade(considerando as tomografias realizadas no trauma)
        #volume por pixel = 0,298 mm³
        
        #funcoes auxiliares
        #Verifica se tem uma região adjacente
        def regiaoAdjacenteY(x,y,z):
            ra = self.regiaoPorPixel.get((z,y-1,x), -1);
            if ra != -1:
                return ra;
            return -1;
    
        def regiaoAdjacenteZ(x,y,z):
            ra = self.regiaoPorPixel.get((z-1,y,x), -1);
            if ra != -1:
                return ra;
            return -1;
        
        def testarCriterio(p):
            minimo = MINIMA_DENSIDADE_SANGUE + 1024;
            maximo = MAXIMA_DENSIDADE_SANGUE + 1024;    
            if p >= minimo and p <= maximo:    
                return True;
            return False;    
        
        def testarOsso(p):
            minimoUHOsso = 95 + 1024;
            if p >= minimoUHOsso:
                return True;
            return False;        
    
        def testarProximidadeOsso(corte, x, y):
            try:
                if testarOsso(corte[y-1][x]):
                    return True;
                if testarOsso(corte[y-1][x-1]):
                    return True;
                if testarOsso(corte[y][x-1]):
                    return True;
                if testarOsso(corte[y+1][x-1]):
                    return True;
                if testarOsso(corte[y+1][x]):
                    return True;
                if testarOsso(corte[y+1][x+1]):
                    return True;
                if testarOsso(corte[y][x+1]):
                    return True;
                if testarOsso(corte[y-1][x+1]):
                    return True;
            except :                
                return True;
            return False;
    
        def unirRegioes(r1, r2):
            #Adiciona todos os pixels de r2 a r1 e apaga r2
            #Antes disso, eu puxo todas as regioes que sao contiguas a r1
            if r1 == r2 or r1 == -1 or r2 == -1: 
                return;
            for p2 in self.regioes[r2].pontos:
                self.regiaoPorPixel[p2] = r1;
                self.regioes[r1].pontos.append(p2);
            del self.regioes[r2];
    
        regiaoAtual = -1;
        indiceRegiao = 0;
        
        t0 = time.perf_counter(); print ('Etapa 1 - Iniciando a Selecao de regioes');
        for z in range (0, len(self.cortesModificados)):
            corteModificado = self.cortesModificados[z];     
            
            #Tira do corte as linhas e colunas que são muito finas, <= SEQUENCIA_MINIMA pixels
            #as quais provavelmente não influenciarão no resultado
            #retirando as linhas curtas
            
            contadorLinha = 0;
            for y in range(0, self.altura):
                for x in range(0, self.largura):
                    pixelHU = corteModificado[y][x];
                    #testa se o pixel atende ao criterio desejado(sangue)          
                    if testarCriterio(pixelHU):
                        #se atende, eu continuo contando
                        contadorLinha += 1;
                    else:
                        #se não atende, eu testo pra ver se é maior do que a sequencia minima
                        if contadorLinha <= SEQUENCIA_MINIMA and contadorLinha > 0:
                            #se não for, eu apago os pixels
                            for apagaPixel in range(x-contadorLinha, x):
                                corteModificado[y][apagaPixel]= MINIMA_DENSIDADE_SANGUE-1+1024;
                        contadorLinha = 0;        
    
            contadorColuna = 0;
            for x in range(0, self.largura):
                for y in range(0, self.altura):
                    pixelHU = corteModificado[y][x];
                    #testa se o pixel atende ao criterio desejado(sangue)          
                    if testarCriterio(pixelHU):
                        #se atende, eu continuo contando
                        contadorColuna += 1;
                    else:
                        #se não atende, eu testo pra ver se é maior do que a sequencia minima
                        if contadorColuna <= SEQUENCIA_MINIMA and contadorColuna > 0:
                            #se não for, eu apago os pixels
                            for apagaPixel in range(y-contadorColuna, y):
                                corteModificado[apagaPixel][x]= MINIMA_DENSIDADE_SANGUE-1+1024;
                        contadorColuna = 0;
            
            #Agora faz a seleção propriamente dita
            for y in range(0, self.altura):
                for x in range(0, self.largura):
                    #Verifica se o pixel atende o criterio de alguma regiao
                    pixelHU = corteModificado[y][x];
                    #testa se o pixel atende ao criterio desejado(sangue)          
                    if testarCriterio(pixelHU) and not testarProximidadeOsso(self.cortesOriginais[z], x, y):
                        if regiaoAtual == -1:
                            #Primeiro pixel dessa regiao
                            indiceRegiao += 1;
                            regiaoAtual = indiceRegiao;
                            self.regioes[regiaoAtual] = Regiao(regiaoAtual, TipoRegiao.PSEUDO_SANGUE, len(self.cortesOriginais), self.cortesOriginais, self.primeiroCorte);
                        #adiciona o ponto a regiao
                        self.regiaoPorPixel[(z,y,x)] = regiaoAtual;
                        self.regioes[regiaoAtual].pontos.append((z,y,x));

                        rAdjacenteY = regiaoAdjacenteY(x, y, z);
                        rAdjacenteZ = regiaoAdjacenteZ(x, y, z);

                        if rAdjacenteY != -1 and rAdjacenteY != regiaoAtual:
                            unirRegioes(rAdjacenteY, regiaoAtual);
                            regiaoAtual = rAdjacenteY;
                        elif rAdjacenteZ != -1 and rAdjacenteZ != regiaoAtual:
                            unirRegioes(rAdjacenteZ, regiaoAtual);
                            regiaoAtual = rAdjacenteZ;                            
                    else:
                        regiaoAtual = -1;
            #t11 = time.perf_counter() - t01; print 'Etapa 1.1 - Análise dos pixels do corte %d: %f segundos' % (z+self.primeiroCorte, t11);      
        t1 = time.perf_counter() - t0; print ('Etapa 1 - Selecao de regioes - finalizada: %f segundos' % (t1));
        
    def processarRegioes(self, criterios):
        t1 = time.perf_counter();

        possuiSangue = False;     
        
        for r in self.regioes:
            regiao = self.regioes[r];
            #Primeiro prepara a regiao para processamento
            regiao.prepararProcessamento();
            resultadoGeralRegiao = True;
            for criterio in criterios:
                resultadoCriterio = criterio.testarCriterio(self, regiao);
                regiao.adicionarResultado(resultadoCriterio);
                if not resultadoCriterio.passou:
                    resultadoGeralRegiao = False;
                    break;
                
            if resultadoGeralRegiao:
                regiao.tipo = TipoRegiao.SANGUE;
                possuiSangue = True;
                
        t2 = time.perf_counter() - t1; print('Etapa 2 - Processamento dos critérios geométricos - finalizada: %f segundos' % (t2));
        
        if possuiSangue:
            #A tomografia foi classificada como AVCh
            self.setDiagnostico('AVCh');
        else:
            self.setDiagnostico('Sem Alteracoes'); 
        
        return possuiSangue;
        
    def getDiagnostico(self):
        return self.diagnostico;
    
    def setDiagnostico(self, d):
        self.diagnostico = d;
class Regiao:
        
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
        
def processarRegiaoTomografia(r, tomografia, regioes, criterios):
    regiao = regioes.get(r);
    #Primeiro prepara a regiao para processamento
    regiao.prepararProcessamento();
    resultadoGeralRegiao = True;
    possuiSangue = False;
    for criterio in criterios:
        resultadoCriterio = criterio.testarCriterio(tomografia, regiao);
        regiao.adicionarResultado(resultadoCriterio);
        if not resultadoCriterio.passou:
            resultadoGeralRegiao = False;
            break;
        
    if resultadoGeralRegiao:
        regiao.tipo = TipoRegiao.SANGUE;
        possuiSangue = True;   

    if possuiSangue:
        #A tomografia foi classificada como AVCh
        tomografia.setDiagnostico('AVCh');
    else:
        tomografia.setDiagnostico('Sem Alteracoes'); 

def processarRegioesTomografia(tomografia, regioes, criterios):
    # numeroProcessadores = 8;
    # chunksize = int(len(regioes)/numeroProcessadores);
    # pool = multiprocessing.Pool(numeroProcessadores);

    t0 = time.perf_counter();

    print('Processando %d regioes' % (len(regioes)));

    for r in regioes:
        processarRegiaoTomografia(r, tomografia, regioes, criterios);

    # ps = [];
    # for r in regioes:
    #     p = multiprocessing.Process(target=processarRegiaoTomografia, args=(r, tomografia, regioes, criterios));
    #     p.start();
    #     ps.append(p);

    # for p in ps:
    #     p.join();

    # pool.imap_unordered(partial(processarRegiaoTomografia, tomografia=tomografia, 
    #                     regioes=regioes, criterios=criterios),
    #                     list(regioes.keys()), chunksize); 
    
    #for r in regioes.keys():
    #    processarRegiaoTomografia(r, tomografia, regioes, criterios);

    # pool.imap(partial(processarRegiaoTomografia, tomografia=tomografia, 
    #                     regioes=regioes, criterios=criterios),
    #                     list(regioes.keys()), chunksize); 

    # pool.close();
    # pool.join();
    
    print('Etapa 2 - Processamento dos critérios geométricos - finalizada: %f segundos' % (time.perf_counter()-t0));

def salvarCorteProcessadoPng(c, tomografia, janelaMinimo, janelaMaximo, windowWidth, 
                             altura, largura, novoCaminho, nomePaciente, primeiroCorte):
    corte = tomografia.getCorteOriginal(c);
    regiaoPorPixel = tomografia.getRegiaoPorPixel();
    regioes = tomografia.getRegioes();
    image_2d = np.array(corte).astype(float);    
    image_2d_scaled = ((np.minimum(np.maximum(image_2d, janelaMinimo), janelaMaximo)  - janelaMinimo) / windowWidth) * 255.0;    
    image_2d_scaled = np.uint8(np.round(image_2d_scaled+0.01));
    pixels = [];
    for y in range(0, altura):
        linha = [];
        for x in range(0, largura):
            #print "Processando corte, x, y (%d, %d, %d)" % (c, x, y);
            r = regiaoPorPixel.get((c,y,x), -1); 
            if  r != -1:
                regiao = regioes[r];
                if regiao.tipo == TipoRegiao.SANGUE:
                    linha.extend([image_2d_scaled[y][x], 30, 30]);
                if regiao.tipo == TipoRegiao.PSEUDO_SANGUE:
                    linha.extend([0, image_2d_scaled[y][x], 0]);
                if regiao.tipo == TipoRegiao.INDETERMINADA:
                    linha.extend([image_2d_scaled[y][x], image_2d_scaled[y][x], image_2d_scaled[y][x]]);                            
            else:
                linha.extend([image_2d_scaled[y][x], image_2d_scaled[y][x], image_2d_scaled[y][x]]);
        pixels.append(linha);    
    arquivoImagem = novoCaminho+'/%s corte processado %03d.png' % (nomePaciente, c+primeiroCorte+1);
    if os.path.exists(arquivoImagem):
        os.remove(arquivoImagem);
    with open(arquivoImagem, 'wb') as png_file:
        w = png.Writer(largura, altura, greyscale=False);
        w.write(png_file, pixels);            

def salvarCorteLimpoPng(c, tomografia, janelaMinimo, janelaMaximo, windowWidth, 
                             altura, largura, novoCaminho, nomePaciente, primeiroCorte
                             ):                             
    corte = tomografia.getCorteOriginal(c);
    image_2d = np.array(corte).astype(float);    
    image_2d_scaled = ((np.minimum(np.maximum(image_2d, janelaMinimo), janelaMaximo)  - janelaMinimo) / windowWidth) * 255.0;    
    image_2d_scaled = np.uint8(np.round(image_2d_scaled+0.01));
    pixels = [];   
    for y in range(0, altura):
        linha = [];
        for x in range(0, largura):
            #print "Processando corte, x, y (%d, %d, %d)" % (c, x, y);
            linha.extend([image_2d_scaled[y][x], image_2d_scaled[y][x], image_2d_scaled[y][x]]);
        pixels.append(linha);    
    arquivoImagem = novoCaminho+'/%s corte limpo %03d.png' % (nomePaciente,
         c+primeiroCorte+1);
    if os.path.exists(arquivoImagem):
        os.remove(arquivoImagem);
    with open(arquivoImagem, 'wb') as png_file:
        w = png.Writer(largura, altura, greyscale=False);
        w.write(png_file, pixels);      
    
def gerarRelatorio(tomografia, regioes):
    pastaPaciente = tomografia.getNomePaciente()+' '+tomografia.getPatientID()+' '+\
                    tomografia.getStudyId()+' ('+tomografia.getStudyDate()+' '+\
                    tomografia.getStudyTime()+') (SATSCAN=%s)'% tomografia.getDiagnostico();
    tomografia.setPaciente(pastaPaciente);
    novoCaminho = 'c:/tomografias/resultados/'+tomografia.getPastaPaciente();
    if not os.path.exists(novoCaminho):
        os.makedirs(novoCaminho);     

    #Imprimo as informacoes gerais da tomografia junto com os cortes sem alteração
    t1 = time.perf_counter();
    print('Etapa 3 - Geração do relatório do processamento da tomografia');
    #Salva as imagens sem processamento e processadas
    # pTomografiaLimpa = multiprocessing.Process(target=salvarTomografiaLimpaPng, args=(tomografia,));
    # pTomografiaLimpa.start();
    # pTomografiaLimpa.join();

    salvarTomografiaLimpaPng(tomografia);

    # pTomografiaProcessada = multiprocessing.Process(target=salvarTomografiaProcessadaPng, args=(tomografia,));
    # pTomografiaProcessada.start();
    # pTomografiaProcessada.join();

    salvarTomografiaProcessadaPng(tomografia);

    #salva o Relatório em HTML
    gerarRelatorioHTML(tomografia, regioes);
    print('Etapa 3 - Geração do relatório do processamento da tomografia - finalizado: %f segundos' % (time.perf_counter() - t1));

def gerarRelatorioHTML(tomografia, regioes):
    #abre o arquivo de relatorio e insere os cabecalhos html
    cInicial = 0;
    qtd = len(tomografia.getCortesOriginais())-1;
    qtdCortes = len(tomografia.getCortesOriginais());
    if qtdCortes - (cInicial + qtd) < qtd and qtdCortes - (cInicial + qtd) != 0:
        qtd = qtdCortes - 1 - cInicial;
    novoCaminho = 'c:/tomografias/resultados/'+tomografia.getPastaPaciente();
    arquivoResultadoHTML = open(novoCaminho+'/index.html', 'w');
    # textoResultado = "<span style='color:green;'>Verdadeiro</span>" if tomografia.getLaudoRadiologista() == tomografia.getDiagnostico() else "<span style='color:red;'>Falso</span>";
    resultado = f"""
        <html>
        <body'>
            <h1>Processamento da Tomografia de Crânio de : {tomografia.getNomePaciente()}</h1>
            <h3>Laudo do Radiologista presente no local: {tomografia.getLaudoRadiologista()}</h3>
            <h3>Laudo do SATScan: {tomografia.getDiagnostico()}</h3>
            <h4>Resultado: {"<span style='color:green;'>Verdadeiro</span>" if tomografia.getLaudoRadiologista() == tomografia.getDiagnostico() else "<span style='color:red;'>Falso</span>"}</h4>
            <table id='cortes'>
                <tr style='vertical-align: top'>
                    <td id='cortesLimpos' style='cursor: pointer'>""";
    for c in range(cInicial, cInicial+qtd): 
        resultado += "<img class='corte' id='%s' src='%s' style='display: none' onmousemove='atualizarCoordenadas(event)' onmouseout='limparCoordenadas(event)'></img>" % ('limpo%d' % (c+tomografia.getPrimeiroCorte()+1) , '%s corte limpo %03d.png' % (tomografia.getNomePaciente(), c+tomografia.getPrimeiroCorte()+1));   
    resultado += """
                    </td>
                    <td id='cortesProcessados' style='cursor: pointer'>""";
    for c in range(cInicial, cInicial+qtd): 
        resultado += "<img class='corte' id='%s' src='%s' style='display: none' onmousemove='atualizarCoordenadas(event)' onmouseout='limparCoordenadas(event)'></img>" % ('processado%d' % (c+tomografia.getPrimeiroCorte()+1), '%s corte processado %03d.png' % (tomografia.getNomePaciente(), c+tomografia.getPrimeiroCorte()+1));   
    resultado += f"""
                    </td>
                    <td style='vertical-align: top; '>   
                        <div>
                            <span id='cx'>X : </span>&nbsp;<span id='cy'>Y : </span>&nbsp;<span id='cz'>Z : </span>&nbsp;<span id='hu'>HU : </span>
                        </div><br/>
                        <div id='regiao'>
                        </div>                    
                    </td>
                </tr>
            </table>
            <script>
                var x = 0;var y = 0; var z = 0;
                var c = {int(tomografia.getPrimeiroCorte()+1)};
                var primeiroCorte = c;
                var cInicial = {int(tomografia.getPrimeiroCorte()+1)};
                var cFinal = {int(tomografia.getPrimeiroCorte()+qtd)};

                function atualizarCorte() {{
                    document.getElementById('limpo'+c).setAttribute('style', 'display: block');
                    document.getElementById('processado'+c).setAttribute('style', 'display: block');                
                }}

                function esconderImagem(c) {{
                    document.getElementById('limpo'+c).setAttribute('style', 'display: none');
                    document.getElementById('processado'+c).setAttribute('style', 'display: none');                
                }}

                function atualizarCoordenadas(e) {{
                    x = e.offsetX;
                    y = e.offsetY;
                    z = c;
                    atualizarInformacoes();
                }}

                function atualizarInformacoes() {{
                    document.getElementById('cx').innerHTML = 'X :'+String(x);
                    document.getElementById('cy').innerHTML = 'Y :'+String(y);
                    document.getElementById('cz').innerHTML = 'Z :'+String(c);
                    hu = pontoHU[c-primeiroCorte][y][x]-1024;
                    document.getElementById('hu').innerHTML = 'HU :'+String(hu);
                    r = regiaoPorPixel[String(c-primeiroCorte)+','+String(y)+','+String(x)];
                    textoRegiao = '';
                    if (r) {{
                        textoRegiao = regioes[r];
                    }} 
                    document.getElementById('regiao').innerHTML = textoRegiao;
                }}

                function limparCoordenadas(e) {{
                    document.getElementById('cx').innerHTML = 'X : ';
                    document.getElementById('cy').innerHTML = 'Y : ';
                }}

                function mudarCorte(event) {{
                    esconderImagem(c);
                    if (event.deltaY < 0)
                        c -= (c <= cInicial) ? 0 : 1;
                    else
                        c += (c >= cFinal) ? 0 : 1;
                    atualizarCorte();
                    z = c;
                    atualizarInformacoes();
                }}

            document.getElementById('cortes').addEventListener('wheel', mudarCorte);

            window.onload = function () {{
                atualizarCorte();
            }};
    """;
    #colocando a array de regiao por pixel
    resultado += "            var regiaoPorPixel = {";
    for pixel in tomografia.getRegiaoPorPixel().keys():
        resultado += "\"%s\" : %d, " % (','.join(map(str, pixel)), tomografia.getRegiaoPorPixel().get(pixel));
    resultado += "};\r\n";
    #colocando a array de regioes
    resultado += "            var regioes = {";
    for r in regioes.keys():
        resultado += "%d : \"%s\", " % (r, regioes.get(r).getDescricaoComResultados());
    resultado += "};\r\n";
    #colocando a array de HU de cada ponto
    resultado += f"            var pontoHU = {json.dumps(tomografia.cortesOriginais) }; ";
    resultado += """    
                </script>
        </body>
        </html>
    """;
    arquivoResultadoHTML.write(resultado);
    arquivoResultadoHTML.close();

def salvarTomografiaLimpaPng(tomografia):
    #cria uma imagem para cada corte da regiao
    cInicial = 0;
    qtd = len(tomografia.getCortesOriginais())-1;

    windowCenter, windowWidth = 35, 85;
    janelaMinimo = windowCenter +1024 - windowWidth/2.0;
    janelaMaximo = windowCenter +1024 + (windowWidth-1)/2.0;
    novoCaminho = 'c:/tomografias/resultados/'+tomografia.getPastaPaciente();

    qtdCortes = len(tomografia.getCortesOriginais());
    if qtdCortes < cInicial + qtd:
        return;
    if qtdCortes - (cInicial + qtd) < qtd and qtdCortes - (cInicial + qtd) != 0:
        qtd = qtdCortes-1 - cInicial; 

    print("Corte inicial: %d, corte final: %d " % (cInicial, cInicial+qtd));
    t0 = time.perf_counter();

    # chunksize = 30;
    # numeroProcessadores = 1;
    # pool = multiprocessing.Pool(numeroProcessadores);
    # pool.imap_unordered(partial(salvarCorteLimpoPng, tomografia=tomografia, janelaMinimo=janelaMinimo,
    #     janelaMaximo=janelaMaximo, windowWidth=windowWidth, altura=tomografia.getAltura(),
    #     largura=tomografia.getLargura(), novoCaminho=novoCaminho,
    #     nomePaciente=tomografia.getNomePaciente(), primeiroCorte=tomografia.getPrimeiroCorte()),
    #     range(cInicial, cInicial+qtd), chunksize);
    # pool.close();
    # pool.join();

    for c in range(cInicial, cInicial+qtd):
        salvarCorteLimpoPng(c, tomografia=tomografia, janelaMinimo=janelaMinimo,
         janelaMaximo=janelaMaximo, windowWidth=windowWidth, altura=tomografia.getAltura(),
         largura=tomografia.getLargura(), novoCaminho=novoCaminho,
         nomePaciente=tomografia.getNomePaciente(), primeiroCorte=tomografia.getPrimeiroCorte());

    print('Salvar imagens limpas - finalizada: %.2f segundos' % (time.perf_counter()-t0));     
    
def salvarTomografiaProcessadaPng(tomografia):
    #cria uma imagem para cada corte da regiao
    cInicial = 0;
    qtd = len(tomografia.getCortesOriginais())-1;
    windowCenter, windowWidth = 35, 85;
    janelaMinimo = windowCenter +1024 - windowWidth/2.0;
    janelaMaximo = windowCenter +1024 + (windowWidth-1)/2.0;
    novoCaminho = 'c:/tomografias/resultados/'+tomografia.getPastaPaciente();

    qtdCortes = len(tomografia.getCortesOriginais());
    if qtdCortes < cInicial + qtd:
        return;
    if qtdCortes - (cInicial + qtd) < qtd and qtdCortes - (cInicial + qtd) != 0:
        qtd = qtdCortes-1 - cInicial; 

    print("Corte inicial: %d, corte final: %d " % (cInicial, cInicial+qtd));
    t0 = time.perf_counter();

    # chunksize = 30;
    # numeroProcessadores = 1;
    # pool = multiprocessing.Pool(numeroProcessadores);    
    # pool.imap_unordered(partial(salvarCorteProcessadoPng, tomografia=tomografia, janelaMinimo=janelaMinimo,
    #     janelaMaximo=janelaMaximo, windowWidth=windowWidth, altura=tomografia.getAltura(),
    #     largura=tomografia.getLargura(), novoCaminho=novoCaminho,
    #     nomePaciente=tomografia.getNomePaciente(), primeiroCorte=tomografia.getPrimeiroCorte()),
    #     range(cInicial, cInicial+qtd), chunksize);            
    # pool.close();
    # pool.join();

    for c in range(cInicial, cInicial+qtd):
        salvarCorteProcessadoPng(c, tomografia=tomografia, janelaMinimo=janelaMinimo,
         janelaMaximo=janelaMaximo, windowWidth=windowWidth, altura=tomografia.getAltura(),
         largura=tomografia.getLargura(), novoCaminho=novoCaminho,
         nomePaciente=tomografia.getNomePaciente(), primeiroCorte=tomografia.getPrimeiroCorte());            

    print('Salvar imagens processadas - finalizada: %.2f segundos' % (time.perf_counter()-t0)); 
    
def main(argv=None):

    #Opens the csv file with all the results read from the files
    reports = Reports('H:/rsna-intracranial-hemorrhage-detection/stage_2_train.csv');

    #Abrindo todas as tomografias de uma pasta
    pastaTomografias = 'H:/rsna-intracranial-hemorrhage-detection/stage_2_train';
    pastas = [pastaTomografias+name for name in os.listdir(pastaTomografias) if os.path.isdir(pastaTomografias+'./'+name)];
    #Lendo todos os criterios a serem utilizados na analise
    criterios = [  # CriterioCircularidadePequenasAreas(),
                 CriterioAVChTamanhoRegiao(), CriterioAVChLimiteDimensional(), 
                 CriterioAVChDensidadeMinima(), #CriterioAVChRelacaoAreaContentora(), 
                 CriterioAVChPertencimentoCaixaCraniana(), CriterioAVChEntornoNegativo() ];
    qtdTomografias = 0;
    qtdSemAlteracoes = 0;
    qtdAVCh = 0;    
    qtdAVChCorreto = 0;
    qtdSemAlteracoesCorreto = 0;
    resultadosAnalises = [];

    htmlResultadoGeral = """

    <html>
    <body>
    <h1>Resultado Geral</h1>

    """;

    htmlTomografias = "";

    for pasta in pastas:
        for tc in next(os.walk(pasta))[1]:
            t0 = time.perf_counter();
            #Abrindo a tomografia para ser analisada
            regioesAnalisadas = dict();
            tomografiaAnalisada = Tomografia();
            tomografiaAnalisada.setPasta(pasta+'/'+tc);
            tomografiaAnalisada.setRegioes(regioesAnalisadas);
            try:
                tomografiaAnalisada.inicializar();
            except:
                print("A pasta %s não contém uma tomografia válida" % pasta);
                continue;
            qtdTomografias += 1;
            
            #Defining the report of the CTScan
            laudo = 'Sem Alteracoes';
            if :
                laudo = 'AVCh';            
            tomografiaAnalisada.setLaudoRadiologista(laudo);
            

            #Acrescentando o laudo a contagem das tomografias
            if tomografiaAnalisada.getLaudoRadiologista() == 'AVCh':
                qtdAVCh += 1;
            else:
                qtdSemAlteracoes += 1;
            #selecionando as regioes de interesse, que correspondem a sangue
            tomografiaAnalisada.selecionarRegioesTomografia();
            #processando as regioes para definir o tipo de cada regiao como sangue ou não
            tomografiaAnalisada.processarRegioes(criterios);
            # pRegioes = multiprocessing.Process(target=processarRegioesTomografia, 
            #                         args=(tomografiaAnalisada, regioesAnalisadas, criterios));
            # pRegioes.start();
            # pRegioes.join();
            #ao final do processamento, poderão existir regioes com sangue ou não
            #eu gero um relatório com todas as regiões
            gerarRelatorio(tomografiaAnalisada, regioesAnalisadas);
            acertouDiagnostico = tomografiaAnalisada.getLaudoRadiologista() == tomografiaAnalisada.getDiagnostico();
            veredicto = "Verdadeiro";
            if not acertouDiagnostico:
                veredicto = "Falso";
            resultadoAnalise = '(%s) %s - Laudo Radiologista: %s - Diagnostico SATScan: %s (%f segundos)' % (veredicto, tomografiaAnalisada.nomePaciente, tomografiaAnalisada.laudoRadiologista, tomografiaAnalisada.getDiagnostico(), time.perf_counter()-t0); 
            pastaPaciente = tomografiaAnalisada.getNomePaciente()+' '+tomografiaAnalisada.getPatientID()+' '+\
                            tomografiaAnalisada.getStudyId()+' ('+tomografiaAnalisada.getStudyDate()+' '+\
                            tomografiaAnalisada.getStudyTime()+') (SATSCAN=%s)'% tomografiaAnalisada.getDiagnostico();
            htmlTomografias += "<div class='tomografia "+veredicto+"'><a href='"+pastaPaciente+"/index.html'>"+resultadoAnalise+"</a></div>";
            print(resultadoAnalise);
            print("");
            resultadosAnalises.append(resultadoAnalise);
            if tomografiaAnalisada.laudoRadiologista == tomografiaAnalisada.getDiagnostico():
                if tomografiaAnalisada.getDiagnostico() == 'AVCh':
                    qtdAVChCorreto += 1;
                else:
                    qtdSemAlteracoesCorreto += 1;
            del tomografiaAnalisada;
            del regioesAnalisadas;
            gc.collect();
            print("---------------------------------------------------------------");
            
    sumario = "<div class='sumario'>\n";
    sumario += "*** Fim do processamento ***</br>\n";
    sumario += "  Resultado Final:</br>\n";
    sumario += "     %d Tomografias</br>\n" % qtdTomografias;
    sumario += "      - %d AVCh</br>\n" % (qtdAVCh);
    sumario += "      - %d Sem alterações</br>\n" % qtdSemAlteracoes;
    sumario += "-----------------------</br>\n";
    sumario += "      - %d AVCh corretamente diagnosticado</br>\n" % qtdAVChCorreto;
    sumario += "      - %d Sem alterações corretamente diagnosticado</br>\n" % qtdSemAlteracoesCorreto;
    sumario += "</div>";

    print(sumario);    
    
    if qtdAVCh > 0:
        sensibilidade = 1.0*qtdAVChCorreto/qtdAVCh*100;
    else:
        sensibilidade = 0;
    if qtdSemAlteracoes > 0:
        especificidade = 1.0*qtdSemAlteracoesCorreto/qtdSemAlteracoes*100;
    else:
        especificidade = 0; 
    acuracia = 1.0*(qtdAVChCorreto+qtdSemAlteracoesCorreto)/qtdTomografias*100;
        
    estatisticas = "<div class='estatisticas'>\n";
    estatisticas = "*** Estatísticas ***<br/>\n";
    estatisticas += "    Sensibilidade = (%d/%d) %.2f%%<br/>\n" % (qtdAVChCorreto, qtdAVCh, sensibilidade);
    estatisticas += "    Especificidade = (%d/%d) %.2f%%<br/>\n"% (qtdSemAlteracoesCorreto, qtdSemAlteracoes, especificidade);
    estatisticas += "    Acurácia = (%d/%d) %.2f%%<br/>\n" % (qtdAVChCorreto+qtdSemAlteracoesCorreto, qtdTomografias, acuracia);
    estatisticas += "</div>";

    print(estatisticas);
    for r in resultadosAnalises:
        print("  "+r);  

    resultadoHtml = open('c:/tomografias/resultados/resultado Geral.html', 'w');
    resultadoHtml.write(htmlResultadoGeral);
    resultadoHtml.write(sumario);
    resultadoHtml.write(estatisticas);
    resultadoHtml.write(htmlTomografias);
    resultadoHtml.close();           
    
if __name__ == '__main__':
    main();