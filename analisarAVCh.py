# encoding: utf-8
# pylint: disable=unnecessary-semicolon
'''
satscan.abrirTomografiaDycom -- shortdesc

@author:     Augusto Camarotti

@license:    GPL 3.0

@contact:    augustoccl@gmail.com
'''

import os;
import time;
import png;
import numpy as np;
import pydicom;
import json;
from criteriosAVCh import CriterioAVChTamanhoRegiao;
from criteriosAVCh import CriterioAVChLimiteDimensional;
from criteriosAVCh import CriterioAVChDensidadeMinima;
from criteriosAVCh import CriterioAVChRelacaoAreaContentora;
from criteriosAVCh import CriterioAVChPertencimentoCaixaCraniana;
from criteriosAVCh import CriterioAVChEntornoNegativo;
from tomography import Tomography;
from resultReports import Reports;
from multiprocessing import Manager, Value, Array;
from multiprocessing.managers import BaseManager;
from functools import partial;
import gc;

from multiprocessing import freeze_support
from _ast import Or

def Enum(*sequential, **named):
    Enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), Enums)

diagnosticoTomografia = Enum('SEM_ALTERACOES', 'AVCh');

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

    #Opens the csv file with all the results related to the slices
    reports = Reports('H:/rsna-intracranial-hemorrhage-detection/stage_2_train.csv');

    #Opening all the slices from a folder
    pastaTomografias = 'H:/rsna-intracranial-hemorrhage-detection/stage_2_train';
    #pastas = [pastaTomografias+name for name in os.listdir(pastaTomografias) if os.path.isdir(pastaTomografias+'./'+name)];
    #Lendo todos os criterios a serem utilizados na analise
    criterios = [  # CriterioCircularidadePequenasAreas(),
                 CriterioAVChTamanhoRegiao(), CriterioAVChLimiteDimensional(), 
                 CriterioAVChDensidadeMinima(), #CriterioAVChRelacaoAreaContentora(), 
                 CriterioAVChPertencimentoCaixaCraniana(), CriterioAVChEntornoNegativo() ];
    slices = 0;
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
    iteratorFolder = os.walk(pastaTomografias);
    #process each slice of the folder
    for dcmFilename in next(iteratorFolder)[2]:
        t0 = time.perf_counter();
        #Abrindo a tomografia para ser analisada
        ctscan = Tomography(pastaTomografias+"/"+dcmFilename);
        slices += 1;
        
        #Defining the report of the CTScan
        report = 'Sem Alteracoes';
        sliceID = dcmFilename.partition('.')[0];
        laudo = 'No alterations';
        if sliceID in reports.getReports():
            laudo = 'intracranial blood';            
        ctscan.setLaudoRadiologista(laudo);        

        #Acrescentando o laudo a contagem das tomografias
        if ctscan.getLaudoRadiologista() == 'intracranial blood':
            qtdAVCh += 1;
        else:
            qtdSemAlteracoes += 1;
        #selecionando as regioes de interesse, que correspondem a sangue
        ctscan.selectRegions();
        #processando as regioes para definir o tipo de cada regiao como sangue ou não
        ctscan.processarRegioes(criterios);
        # pRegioes = multiprocessing.Process(target=processarRegioesTomografia, 
        #                         args=(tomografiaAnalisada, regioesAnalisadas, criterios));
        # pRegioes.start();
        # pRegioes.join();
        #ao final do processamento, poderão existir regioes com sangue ou não
        #eu gero um relatório com todas as regiões
        gerarRelatorio(ctscan);
        acertouDiagnostico = ctscan.getLaudoRadiologista() == ctscan.getDiagnostico();
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
    sumario += "     %d Tomografias</br>\n" % slices;
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
    acuracia = 1.0*(qtdAVChCorreto+qtdSemAlteracoesCorreto)/slices*100;
        
    estatisticas = "<div class='estatisticas'>\n";
    estatisticas = "*** Estatísticas ***<br/>\n";
    estatisticas += "    Sensibilidade = (%d/%d) %.2f%%<br/>\n" % (qtdAVChCorreto, qtdAVCh, sensibilidade);
    estatisticas += "    Especificidade = (%d/%d) %.2f%%<br/>\n"% (qtdSemAlteracoesCorreto, qtdSemAlteracoes, especificidade);
    estatisticas += "    Acurácia = (%d/%d) %.2f%%<br/>\n" % (qtdAVChCorreto+qtdSemAlteracoesCorreto, slices, acuracia);
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