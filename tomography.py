import time;
import region;
import copy;
import pydicom;

MINIMA_DENSIDADE_SANGUE = 48;
MAXIMA_DENSIDADE_SANGUE = 78; 
SEQUENCIA_MINIMA = 4;   

class Tomography:

    originalSlice = None;
    regioes = None;
    regiaoPorPixel = None;
    pastaPaciente = None;
    diagnostico = None;
    regiaoPorPixel = {};
   
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

    def __init__(self, dcmFilename):
        self.largura = None;
        self.altura = None;
        self.diagnostico = "";
        ds = pydicom.dcmread(dcmFilename);
        self.largura, self.altura = ds.pixel_array.shape;
        self.nomePaciente = str(ds.PatientID);

        self.slice = ds.pixel_array.tolist();
        self.originalSlice = copy.deepcopy(self.slice); 

        print("Slice of %s succesfully read" % self.nomePaciente);

    #auxilary functions
    def regiaoAdjacenteY(self, x,y):
        ra = self.regiaoPorPixel.get((y-1,x), -1);
        if ra != -1:
            return ra;
        return -1;

    def regiaoAdjacenteZ(self, x,y,z):
        ra = self.regiaoPorPixel.get((y,x), -1);
        if ra != -1:
            return ra;
        return -1;
    
    def testarCriterio(self, p):
        minimo = MINIMA_DENSIDADE_SANGUE + 1024;
        maximo = MAXIMA_DENSIDADE_SANGUE + 1024;    
        if p >= minimo and p <= maximo:    
            return True;
        return False;    
    
    def testarOsso(self, p):
        minimoUHOsso = 95 + 1024;
        if p >= minimoUHOsso:
            return True;
        return False;        

    def testarProximidadeOsso(self, corte, x, y):
        try:
            if self.testarOsso(corte[y-1][x]):
                return True;
            if self.testarOsso(corte[y-1][x-1]):
                return True;
            if self.testarOsso(corte[y][x-1]):
                return True;
            if self.testarOsso(corte[y+1][x-1]):
                return True;
            if self.testarOsso(corte[y+1][x]):
                return True;
            if self.testarOsso(corte[y+1][x+1]):
                return True;
            if self.testarOsso(corte[y][x+1]):
                return True;
            if self.testarOsso(corte[y-1][x+1]):
                return True;
        except :                
            return True;
        return False;

    def unirRegioes(self, r1, r2):
        #Adiciona todos os pixels de r2 a r1 e apaga r2
        #Antes disso, eu puxo todas as regioes que sao contiguas a r1
        if r1 == r2 or r1 == -1 or r2 == -1: 
            return;
        for p2 in self.regioes[r2].pontos:
            self.regiaoPorPixel[p2] = r1;
            self.regioes[r1].pontos.append(p2);
        del self.regioes[r2];        
    
    
    def selectRegions(self):
        #===========================================================================
        # Selecting regions which could possibly mean intracranial blood
        #===========================================================================
        print ("Analyzing tomography slice %s - %s" % (self.nomePaciente, self.laudoRadiologista));
                
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
        
   
        regiaoAtual = -1;
        indiceRegiao = 0;
        
        t0 = time.perf_counter(); print ('Etapa 1 - Iniciando a Selecao de regioes');
        
        #Tira do corte as linhas e colunas que são muito finas, <= SEQUENCIA_MINIMA pixels
        #as quais provavelmente não influenciarão no resultado
        #retirando as linhas curtas        
        contadorLinha = 0;
        for y in range(0, self.altura):
            for x in range(0, self.largura):
                pixelHU = self.slice[y][x];
                #testa se o pixel atende ao criterio desejado(sangue)          
                if self.testarCriterio(pixelHU):
                    #se atende, eu continuo contando
                    contadorLinha += 1;
                else:
                    #se não atende, eu testo pra ver se é maior do que a sequencia minima
                    if contadorLinha <= SEQUENCIA_MINIMA and contadorLinha > 0:
                        #se não for, eu apago os pixels
                        for apagaPixel in range(x-contadorLinha, x):
                            self.slice[y][apagaPixel]= MINIMA_DENSIDADE_SANGUE-1+1024;
                    contadorLinha = 0;        

        contadorColuna = 0;
        for x in range(0, self.largura):
            for y in range(0, self.altura):
                pixelHU = self.slice[y][x];
                #testa se o pixel atende ao criterio desejado(sangue)          
                if self.testarCriterio(pixelHU):
                    #se atende, eu continuo contando
                    contadorColuna += 1;
                else:
                    #se não atende, eu testo pra ver se é maior do que a sequencia minima
                    if contadorColuna <= SEQUENCIA_MINIMA and contadorColuna > 0:
                        #se não for, eu apago os pixels
                        for apagaPixel in range(y-contadorColuna, y):
                            self.slice[apagaPixel][x]= MINIMA_DENSIDADE_SANGUE-1+1024;
                    contadorColuna = 0;
        
        #Agora faz a seleção propriamente dita
        for y in range(0, self.altura):
            for x in range(0, self.largura):
                #Verifica se o pixel atende o criterio de alguma regiao
                pixelHU = self.slice[y][x];
                #testa se o pixel atende ao criterio desejado(sangue)          
                if self.testarCriterio(pixelHU) and not self.testarProximidadeOsso(self.originalSlice, x, y):
                    if regiaoAtual == -1:
                        #Primeiro pixel dessa regiao
                        indiceRegiao += 1;
                        regiaoAtual = indiceRegiao;
                        self.regioes[regiaoAtual] = region.Region(regiaoAtual, region.TipoRegiao.PSEUDO_SANGUE, len(self.cortesOriginais), self.cortesOriginais, self.primeiroCorte);
                    #adiciona o ponto a regiao
                    self.regiaoPorPixel[(z,y,x)] = regiaoAtual;
                    self.regioes[regiaoAtual].pontos.append((z,y,x));

                    rAdjacenteY = self.regiaoAdjacenteY(x, y);
                    rAdjacenteZ = self.regiaoAdjacenteZ(x, y);

                    if rAdjacenteY != -1 and rAdjacenteY != regiaoAtual:
                        self.unirRegioes(rAdjacenteY, regiaoAtual);
                        regiaoAtual = rAdjacenteY;
                    elif rAdjacenteZ != -1 and rAdjacenteZ != regiaoAtual:
                        self.unirRegioes(rAdjacenteZ, regiaoAtual);
                        regiaoAtual = rAdjacenteZ;                            
                else:
                    regiaoAtual = -1;

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
                regiao.tipo = region.TipoRegiao.SANGUE;
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
