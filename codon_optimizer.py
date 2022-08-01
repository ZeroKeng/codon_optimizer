from operator import length_hint
from random import seed
import pandas as pd
from scvelo import DataFrame
import streamlit as st
import numpy as np
import seaborn as sns
import base64
import matplotlib.pyplot as plt
from PIL import Image
import time
from dnachisel import *
from Bio import Seq
import re
import random
np.random.seed(123)
random.seed(123)
# 正文图
image = Image.open('images.png')
st.image(image, use_column_width=True)
#

st.markdown("""
            # Codon Optimizer
            
            """)
st.sidebar.header('Input Parameters')
sequence_name = st.sidebar.text_input('Sequence name')
st.header('')


#DNAchisel
species = st.sidebar.selectbox('Species to optimise for',('h_sapiens','m_musculus','s_cerevisiae','e_coli'))
min_gc, max_gc = st.sidebar.slider('GC content',0.0,1.0,(0.4,0.6))
n_cgc = int(st.sidebar.number_input('Avoid continuous G/C number',3)) #number of constitute G/C to avoid
window_l = int(st.sidebar.number_input('Sliding window length',50)) #number of constitute G/C to avoid
kmer = int(st.sidebar.number_input('Avoid continuous kmer length',4))

sequence = st.text_area('Input sequence')
sequence = re.sub(r'[^a-zA-Z]', '', sequence).upper()

start, end = st.sidebar.slider('Translation start site',1,len(sequence),(1,len(sequence)))
K = st.sidebar.number_input('Kmer length to plot kmer matrix',5)
#print(start,end)

def solveProblem():
    problem = DnaOptimizationProblem(
        sequence = sequence,
        constraints=[
            AvoidPattern('{}xC'.format(n_cgc)),
            AvoidPattern('{}xG'.format(n_cgc)),
            AvoidPattern('2x{}mer'.format(kmer)),
            #AvoidPattern(r'([ATCG]{' + str(kmer) + ',})[ATCG]+\1{2,}'),
            EnforceGCContent(mini=min_gc, maxi=max_gc , window=window_l),
            EnforceTranslation(location=(start-1, end))
        ],
        objectives=[CodonOptimize(species=species)]
    )

    problem.resolve_constraints()
    problem.optimize()
    st.write(problem.constraints_text_summary())
    st.write(problem.objectives_text_summary())
    final_sequence = problem.sequence
    st.text_area('Optimized sequence',final_sequence)
    return final_sequence

def filedownload(final_sequence):
        b64 = base64.b64encode(final_sequence.encode()).decode()  # strings <-> bytes conversions
        href = f'<a href="data:file/fasta;base64,{b64}" download="{sequence_name}.fasta">Download FASTA File</a>'
        return href
    
def analysis_progress():
    my_bar = st.progress(0)
    for percent_complete in range(100):
        time.sleep(0.01)
        my_bar.progress(percent_complete + 1)
  
def plotGC(seq,window_l):
    st.header('GC content in each window')
    GC_list = []
    for i in range(len(seq) - window_l):
        GC = len(re.findall(r'[GC]',seq[i:i + window_l + 1]))/window_l
        GC_list.append(GC)
    st.line_chart(pd.DataFrame({'GC content':GC_list}))

def plotKmerMatrix(seq,K):
    kmer_list = [seq[i:i+K] for i in range(len(seq)-K +1)]
    kmer_matrix = np.zeros((len(kmer_list),len(kmer_list)))
    for i in range(len(kmer_list)):
        for j in range(len(kmer_list)):
            if kmer_list[i] == kmer_list[j]:
                kmer_matrix[i,j] = 1
            else:
                kmer_matrix[i,j] = 0
    with sns.axes_style():
        f, ax = plt.subplots(figsize=(7, 5))
        ax = sns.heatmap(kmer_matrix,cbar=False)
    st.pyplot(f)
        

if st.button('Start optimizing'):
    analysis_progress()
    try:
        final_sequence = solveProblem()
        #print(final_sequence)  
        st.markdown(filedownload(final_sequence), unsafe_allow_html=True)
        st.success('Optimization is successful!')
        plotGC(sequence,window_l)
        plotGC(final_sequence,window_l)
        st.header('Kmer correlation matrix heatmap')
        st.write('Kmer correlation matrix heatmap shows the local identity across the whole optimized sequence')
        st.write('Kmer correlation matrix heatmap before optimization')
        plotKmerMatrix(sequence,K)
        st.write('Kmer correlation matrix heatmap after optimization')
        plotKmerMatrix(final_sequence,K)
        
    except(NoSolutionError):
        st.write('This sequence cannot be optimized under current parameters')
        st.error('Optimization is failed!')
        
    except(ValueError):
         st.error('Location length should be multiple of three!')

    