import streamlit as st
import pandas as pd
from process_safety_utils import color_picker,main, summary,create_dataframe
from config import hefg_list
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

st.markdown("# Process Safety")
st.markdown('**Note - Please do not post target or intermediate structure information externally**.')



smile = st.text_input('Enter Smile Code', value= 'O1N=C(C=C1)C(C1)=CC(=CC=1N=[N+]=[N-])CC')
# formula = st.text_input('Enter Molecule Formula', 'C11H10N4O')




df = pd.read_excel('atomic masses.xlsx')

mass_df = df[df.columns[[1, 3]]]

sym_to_mass = {}

for i in range(mass_df.shape[0]):
    sym = mass_df.iloc[i, 0]
    mass = mass_df.iloc[i, 1]
    sym_to_mass[sym] = mass




# st.markdown("""
#     <style>
#     .big-font {
#         font-size:19px !important;
#         font-family: serif !important;
#     }
#     </style>
#     """, unsafe_allow_html=True)

m = st.markdown("""
<style>
div.stButton > button:first-child {
    background-color: #0C1B2A;
    color:#ffffff;
    border:None;
}
div.stButton > button:first-child:focus {
    background-color: #0C1B2A;
    color:#ffffff;
    border:None;
}
</style>""", unsafe_allow_html=True)




if st.button('Calculate'):
    with st.spinner('Processing...'):
        mol = Chem.MolFromSmiles(str(smile))
        formula = CalcMolFormula(mol)
        hefg, r_v, oxy, group = main(smile, hefg_list, sym_to_mass, formula)


        color, text = color_picker(oxy)
        df = create_dataframe(group, r_v, oxy, text)
        st.download_button(
            label= 'Download Data',
            data = df,
            file_name = 'project_safety.csv',
            mime= 'text/csv'
        )
        col1, col2, col3, col4 = st.columns((1.3, 1, 1, 1))
        with col1:
            st.text('HEFG')
            for i in range(len(group)):
                if i + 1 < len(group):
                    st.subheader(str(group[i] + ','))
                else:
                    st.subheader(str(group[i]))
        
        with col2:
            st.text('Rule Six')
            st.subheader(round(r_v, 2))
        with col3:
            st.text('Oxygen Balance')
            st.subheader(oxy)
        with col4:
            st.text('Hazard Rank')
            st.subheader(text)
        
        st.text('')
        st.text('Summary')
        summary(group, r_v, oxy, text) 


# s = ['C1OC1', 'C1CNC1', 'OO', 'N=[N+]=[N-]', 'N1=CC=CO1']
# # smile, formula = 'O1N=C(C=C1)C(C1)=CC(=CC=1N=[N+]=[N-])CC' , 'C11H10N4O'
# # smile, formula = 'C1N(CC1)CC1=CC(=CC=C1)Br', 'C10H12BrN'
# smile, formula = 'O1C(C1)C1=CC(=CC=C1)C(OO)=O', 'C9H8O4'




# add summary below to display why hazard rank is high or low
