def fasta_downloader(id_coati):
    "fasta_downloader permite imgresar comandos para busqueda en la web y descargarlo al dispositivo conectado"  
    from Bio import Entrez
    in_sequence = open(id_coati, "r")
    out_sequence_gb = open("data/coati.gb", "w")
    out_sequence_fasta = open("data/coati.fasta", "w")
    
    for linea in in_sequence:
        Entrez.email = "omar.caiza@est.ikiam.edu.ec"
        handle = Entrez.efetch(db = "nucleotide", id = linea, rettype = "fasta", retmode = "text")
        data = (handle.read())
        out_sequence_fasta.write(data)
    out_sequence_fasta.close()
    
    for linea in in_sequence:
        Entrez.email = "omar.caiza@est.ikiam.edu.ec"
        handle = Entrez.efetch(db = "nucleotide", id = linea, rettype = "gb", retmode = "text")
        data = (handle.read())
        out_sequence_gb.write(data)
    out_sequence_gb.close()
#----------------------------------------------------------------------------------------------------#
def aligment (archivo_fasta):
    "reconocimiento y alineamiento de las secuencias desargadas y guardadas en las carpetas coati.dnd y coati.aln"
    from Bio.Align.Applications import ClustalwCommandline
    import os
    from Bio import AlignIO
    from Bio import Phylo
    clustalw_exe = r"C:\Program files (x86)\ClustalW2\clustalw2.exe"
    clustalw_cline = ClustalwCommandline(clustalw_exe, infile = "data\coati.fasta")
    assert os.path.isfile(clustalw_exe), "Clustal_w executable is missing or not found"
    stdout, stderr = clustalw_cline()
    print (clustalw_cline)
    ClustalAlign = AlignIO.read("data/coati.aln", "clustal")
    print(ClustalAlign)
    tree = Phylo.read("data/coati.dnd", "newick")

#-----------------------------------------------------------------------------------------------------#
def tree (alineacion):
    "la funcion 'tree' grafica el alineamiento"
    from Bio import Phylo 
    from Bio import SeqIO
    from Bio import AlignIO
    from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
    from Bio.Phylo.TreeConstruction import DistanceCalculator
    
    with open (alineacion, "r") as aln:
        alignment = AlignIO.read(aln, "clustal")
        calculator = DistanceCalculator("blosum62")
        distance_matrix = calculator.get_distance(alignment)
        constructor = DistanceTreeConstructor(calculator)
    ##ARBOL
    align_total = constructor.build_tree(alignment)
    align_total.rooted = True
    Phylo.write(align_total, "coati.xml", "phyloxml")
    ## librerias para construir el arbol
    import matplotlib
    import matplotlib.pyplot as plt
    ### configuracion de medidas, colores, etc del arbol
    fig = plt.figure(figsize = (30, 40), dpi = 100)
    matplotlib.rc("font", size = 20)
    matplotlib.rc("xtick", labelsize = 20)           
    matplotlib.rc("ytick", labelsize = 20)
    axes = fig.add_subplot(1, 1, 1)
    Phylo.draw(align_total, axes = axes)
                
    #creacion del archivo coati_phylotree.pdf con el grafico
    fig.savefig("data/coati_phylotree.pdf", dpi = 500 )
