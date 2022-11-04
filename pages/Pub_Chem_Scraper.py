from lib2to3.pgen2 import driver
import streamlit as st
import pandas as pd
from process_safety_utils import *
from pubchem_utils import get_driver, find_cas_number_link, extract_data, make_df, merge_dataframe, create_category_df

st.markdown("# PubChem Scraper")
st.markdown('**Note - Please do not post target or intermediate structure information externally**.')

category_file = 'Categories.xlsx'

category = {'Green' : ['H302','H312','H332','H333','H303','H305', 'H313','H315','H316','H319','H320','H335'],
           'Amber' : ['H301','H311','H331','H300','H310', 'H304', 'H314','H336'],
            'Red' : ['H330','H340','H350','H360', 'H317', 'H334', 'H318', 'H341', 'H351', 'H361', 'H370', 'H371','H372','H373']}





cas_no = st.text_input('CAS Number', '57-27-2,102-54-5, 108-88-3')


cas_no_list = cas_no.split(',')

def data_scraper(start_link, driver):
    link = find_cas_number_link(start_link, driver)
    driver.get(link)
    data = extract_data(category, driver)
    df = make_df(data)
    return df

def category_scraper(data_df, category_file):
    cat_df = create_category_df(data_df, category_file)
    return cat_df


def main():
    dataframe = pd.DataFrame()
    for i in cas_no_list:
        start_link = "https://pubchem.ncbi.nlm.nih.gov/#query=" + str(i)
        

        driver.get(start_link)
        driver.implicitly_wait(2)
        
        data_df = data_scraper(start_link, driver)
        dataframe = merge_dataframe(dataframe, data_df)

    return dataframe


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

if st.button('Search'):
    with st.spinner('Processing...'):
        driver = get_driver() 

        dataframe = main()

        cat_df = category_scraper(dataframe, category_file)

        st.caption('Data')
        st.write(dataframe)
        st.caption('Category')
        st.write(cat_df)

        