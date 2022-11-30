import re
from selenium import webdriver
from selenium.webdriver import ChromeOptions
from selenium.webdriver.chrome.service import Service
from webdriver_manager.chrome import ChromeDriverManager
from selenium.webdriver.common.by import By
import pandas as pd
import streamlit as st


def get_driver():
    option = webdriver.ChromeOptions()
    option.add_argument('--headless')
    option.add_argument('--no-sandbox')
    service = Service(ChromeDriverManager().install())
    driver = webdriver.Chrome(options = option, service= service)
    return driver

def find_cas_number_link(start_link, driver):

    tmp = driver.find_elements(By.CSS_SELECTOR, 'span.breakword')
    driver.implicitly_wait(2)
    cid = tmp[1].text
    n_link = start_link.split('#')[0]
    link = n_link + 'compound/' + cid
    n_link = start_link.split('#')[0]
    link = n_link + 'compound/' + cid
    return link


def get_name(data, driver):
    """To get the name of molecule

    Args:
        data {dict}: empty dict to store name and other information.
        driver {selenium web driver}: to load the website and extract data.

    Returns:
        data {dict}: returns a dictionary of information stored.
    """
    try:
        name = driver.find_element(By.CSS_SELECTOR, 'h1.m-zero')
        compound_name = name.text
        data['Name'] = compound_name
    except:
        # st.write("[INFO] Name is not avialable")
        data['Name'] = 'None'

    return data

def get_summary(data, driver):
    """To get molecule formula and weight

    Args:
        data {dict}: dictionary to store information extracted from website.
        driver {selenium web driver}: to load the website and extract data.

    Returns:
        data {dict}: returns a dictionary of information stored.
    """
    try:
        temp_list = ['Molecular Formula', 'Molecular Weight']
        summary = driver.find_element(By.CSS_SELECTOR, 'table.summary')
        tr_tags = summary.find_elements(By.CSS_SELECTOR, 'tr')
        
        if tr_tags:
            for check_string in temp_list:
                for tags in tr_tags:
                    if check_string in tags.text:
                        s = tags.text.split()
                        data[check_string] = s[2]

                        break
    except:
        # st.write("[INFO] Molecule Formula and Molecule weight not found")
        data['Molecular Formula'] ='None'
        data['Molecular Weight'] = 'None'

    return data

def get_smile(data, driver):
    """To get the smile formula.

    Args:
        data {dict}: dictionary to store information extracted from website.
        driver {selenium web driver}: to load the website and extract data.

    Returns:
        data {dict}: returns a dictionary of information stored.
    """
    try:
        smile = driver.find_element(By.ID, 'Canonical-SMILES')
        smile_code = smile.find_element(By.CSS_SELECTOR, 'p').text
        data['Smile'] = smile_code
    except:
        # st.write('[INFO] Smile not found')
        data['Smile'] = 'None'

    return data

def get_density(data, driver):
    """To get the density.

    Args:
        data {dict}: dictionary to store information extracted from website.
        driver {selenium web driver}: to load the website and extract data.

    Returns:
        data {dict}: returns a dictionary of information stored.
    """
    match_string = '°C'
    try:
        density_tag = driver.find_element(By.ID, 'Density')
        denisty_content = density_tag.find_elements(By.CSS_SELECTOR, 'div.section-content-item')
        for i in range(len(denisty_content)):
            m = denisty_content[i].find_element(By.CSS_SELECTOR, 'p')
            if match_string in m.text:
                data['Density'] = m.text
    except:
        data['Density'] = None

    return data

def get_h_statemenmt(info, hazard):
    """To get hazard statement.

    Args:
        info : information scraped from website
        hazard {dict}: empty dictionary to store hazard statements.
        driver {selenium web driver}: to load the website and extract data.

    Returns:
        hazard {dict}: dictionary to store hazard statements.
    """
    for j in info:
        r = re.findall(r'H[0-9][0-9][0-9].*', str(j.text))
        if r:
            temp = r[0].split(':')
            if '(' in temp[0]:
                z = temp[0].split(' ')
                hazard[z[0]] = temp[1]
            else:
                hazard[temp[0]] = temp[1]
    return hazard

def get_ghs(hazard, driver):
    """To get the information from website of GHS classification.

    Args:
        hazard {dict}: empty dictionary to store hazard statements.
        driver {selenium web driver}: to load the website and extract data.

    Returns:
        hazard {dict}: empty dictionary to store hazard statements.
    """
    try:
        ghs = driver.find_element(By.ID, 'GHS-Classification')
        ghs_temp = ghs.find_elements(By.CSS_SELECTOR, 'div.breakword')
        hazard = get_h_statemenmt(ghs_temp, hazard)

        ghs_string = ghs.find_elements(By.CSS_SELECTOR, 'p')
        hazard = get_h_statemenmt(ghs_string, hazard)
    except:
        pass

    
    return hazard

def make_single(hazard, data):
    """_summary_

    Args:
        hazard {dict}: empty dictionary to store hazard statements.
        data {dict}: dictionary to store information extracted from website.
        driver {selenium web driver}: to load the website and extract data.

    Returns:
        data {dict}: dictionary to store information extracted from website.
    """
    string = ''
    if hazard:
        last_item = list(hazard)[-1]
        for item in hazard:
            if item != last_item:
                if '*' in item:
                    i = item.strip('*')
                    string += str(i) + ' - ' + str(hazard[item]) + str(' , ')
                else:
                    string += str(item) + ' - ' + str(hazard[item]) + str(' , ')
            else:
                if '*' in item:
                    i = item.strip('*')
                    string += str(i) + ' - ' + str(hazard[item])
                else:
                    string += str(item) + ' - ' + str(hazard[item])
        
        data['Hazard'] = string
    else:
        string = 'Not found'
        data['Hazard'] = string
    
    return data

def check_category(category, hazard):
    """_summary_

    Args:
        category (_type_): _description_
        hazard (_type_): _description_
        driver (_type_): _description_

    Returns:
        _type_: _description_
    """
    category_item = {}
    tmp = list(hazard.keys())
    for i in tmp:
        if '(' in i:
            x = i.split(' ')
            y = x[0]
        else:
            y = i
    
        for j in category.keys():
            if y in category[j]:
                if j in category_item.keys():
                    category_item[j].append(y)
                else:
                    category_item[j] = [y]
                
    
    return category_item

def final_category(category, hazard, data):
    """_summary_

    Args:
        category (_type_): _description_
        hazard (_type_): _description_
        data (_type_): _description_
        driver (_type_): _description_

    Returns:
        _type_: _description_
    """
    if hazard:
        category_item = check_category(category, hazard)
        if 'Red' in category_item.keys():
            cat = 'Red'
        elif 'Amber' in category_item.keys():
            cat = 'Amber'
        elif 'Green' in category_item.keys():
            cat = 'Green'
        else:
            cat = 'Special'
        
        data['Category'] = cat
    
    else:
        data['Category'] = 'Not found'
    
    return data

def create_df_data(data):
    """_summary_

    Args:
        data (_type_): _description_

    Returns:
        _type_: _description_
    """

    df = pd.DataFrame({'Name': [data['Name']],
                    'Molecular Formula': [data['Molecular Formula']],
                    'Molecular Weight': [data['Molecular Weight']],
                    'Smile': [data['Smile']],
                    'Density': [data['Density']],
                    'Hazard': [data['Hazard']],
                    'Category': [data['Hazard']]})

    return df

def find_major_cateogory(category_df, category_file):
    """_summary_

    Args:
        category_df (_type_): _description_
        category_file (_type_): _description_

    Returns:
        _type_: _description_
    """
    if 'Red' in category_df:
        cat = pd.read_excel(category_file, sheet_name= 'Red')
    elif 'Amber' in category_df:
        cat = pd.read_excel(category_file, sheet_name= 'Amber')
    elif 'Green' in category_df:
        cat = pd.read_excel(category_file, sheet_name= 'Green')
    else:
        cat = pd.read_excel(category_file, sheet_name= 'Special')
    
    return cat


def extract_data(category, driver):
    data = {}
    hazard = {}
    # driver = get_driver(link)
    data = get_name(data, driver)
    data = get_summary(data, driver)
    data = get_smile(data, driver)
    data = get_density(data, driver)       

    hazard = get_ghs(hazard, driver)
    data = make_single(hazard, data)
    data = final_category(category, hazard, data)

    return data

def make_df(data):
    df = pd.DataFrame([data])
    return df

def merge_dataframe(dataframe, df):

    dataframe = pd.concat([dataframe, df], ignore_index= True)

    return dataframe

def create_category_df(dataframe, category_file):
    category_df = dataframe.loc[:, 'Category']
    cat_df = find_major_cateogory(list(category_df), category_file)

    return cat_df


