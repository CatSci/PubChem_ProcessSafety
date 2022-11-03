import pandas as pd
from collections import defaultdict
import streamlit as st


###############
# To find HEFG
###############

def find_part(data, hefg_list):
    """To find HEFG present and count of HEFG

    Arguments:  
        data {[String]} -- A string of Smile.

        hefg_list {[List]} -- A list of HEFG to check if present in data.

    Returns:
        count {[integer]} -- number of HEFG present

        group {[List]} -- A list of all HEFG found in data.
    """
    count = 0
    group = []
    for i in range(len(hefg_list)):
        if hefg_list[i] in data:
            count += 1
            group.append(hefg_list[i])
    return count, group

####################
# To find 'C' atoms
####################

def rule_six(data, hefg):
    """To calculate rule six

    Args:
        data {[String]} -- A string of Smile.
        hefg {[hefg]} -- Count of HEFG present

    Returns:
        x {[integer]} -- value of rule six 
    """
    count = 0
    for i in range(len(data)):
        if data[i] == 'C':
            count += 1
    if hefg > 0:
        x = count / hefg
    
        if x < 6:
            return x
    else:
        x = 0
    
    return x


#################
# To count atoms
#################

def countOfAtoms(formula: str) -> str:
    """ To calculate number of atoms in a chemical formula string

    Args:
        formula {[String]}: string of chemical formula

    Returns:
        atoms {[List]} -- List of all atoms and the count values of each atoms
    """

    dt = defaultdict(int)
    stack = [1]
    digits = ""
    lowers = ""
    for element in formula[::-1]:
        if element.isdigit():
            digits = element + digits
        elif element.islower():
            lowers = element + lowers
        elif element == ")":
            stack.append(stack[-1]*int(digits or 1))
            digits = ""
        elif element == "(":
            stack.pop()
        #if element is an uppercase letter
        else:
            element = element + lowers
            dt[element] = dt[element]+stack[-1]*(int(digits or 1))
            digits = ""
            lowers = ""
    result = []
    for key, value in sorted(dt.items()):
        if value == 1:
            value = ""
        result.append(key)
        result.append(str(value))
    
    atoms = []
    for i in result:
        if i == '':
            atoms.append('1')
        else:
            atoms.append(i)
    return atoms

###############################
# To calculate molecule weight
###############################

def molecule_weight(sym_to_mass, atoms):
    """ To calculate molecular weight of compound

    Args:
        sym_to_mass {[Dataframe]} -- Dataframe of atoms mass
        atoms {[List]} --  List of all atoms and the count values of each atoms

    Returns:
        mass {[integer]} -- Molecule mass of compound
    """
    mass = 0
    for i in range(len(atoms)):
        if atoms[i] in sym_to_mass:
            m = sym_to_mass[atoms[i]]
            a = int(atoms[i + 1])
            mass =  mass + a*m

    return mass 


####################
# To find oxy params
####################

def oxy_params(atoms):
    """ To calculate params for oxygen balance

    Args:
        atoms {[List]} -- List of all atoms and the count values of each atoms

    Returns:
        x {[Integer]} -- number of Carbon atoms
        y {[Integer]} -- number of Hydrogen atoms
        z {[Integer]} -- number of Oxygen atoms
    """
    for i in range(len(atoms)):
        if 'C' not in atoms:
            x = 0
        if 'H' not in atoms:
            y = 0
        if 'O' not in atoms:
            z = 0

        if atoms[i] == 'C':
            x = int(atoms[i + 1])
        if atoms[i] == 'H':
            y = int(atoms[i + 1])
        if atoms[i] == 'O':
            z = int(atoms[i + 1])
  
    return x, y, z

#############################
# To calculate oxygen balance
#############################

def oxygen_balance(sym_to_mass, formula):
    """ To calculate oxygen balance

    Args:
        sym_to_mass (_type_): _description_
        formula (_type_): _description_


    Returns:
        oxy_bal {[Integer]} -- Value of oxygen balance
    """
    atoms = countOfAtoms(formula)  
    mw = molecule_weight(sym_to_mass, atoms)
    x, y, z = oxy_params(atoms)
    oxy_bal = round(-1600 * (2*x + y/2 -z) / mw)
    return oxy_bal


 

def color_picker(oxy):
    color = ''
    text = ''
    if oxy > 160 or oxy < -240:
        color = '#C2D69B'
        text = 'Low'
    elif 80 < oxy < 160 or -120 > oxy > -240:
        color = '#FFC000'
        text = 'Medium'
    elif 80 > oxy > -120:
        color = '#FF0000'
        text = 'High'

    return color, text

def summary(group, r_v, oxy, text):
    """ To write summary points

    Args:
        group {[List]} -- A list of all HEFG found in smile
        r_v {[Integer]} -- value of rule six
        oxy {[Integer]} -- value of oxygen balance
        text {[String]} -- Hazard rank text
    """
    if 'None' in group:
        st.markdown('<p class="big-font">No High Energy Function groups found.</p>', unsafe_allow_html=True)
    else:
        st.markdown('<p class="big-font">{} High Energy Function groups found.</p>'.format(len(group)), unsafe_allow_html=True)
    if r_v < 6 and r_v != 0:
        st.markdown('<p class="big-font">Failed Rule Six Because its less than 6.</p>'.format(r_v), unsafe_allow_html=True)
    elif r_v >= 6:
        st.markdown('<p class="big-font">Passed Rule Six Because its greater than 6.</p>'.format(r_v), unsafe_allow_html=True)
    else:
        st.markdown('<p class="big-font">Passed Rule Six.</p>'.format(r_v), unsafe_allow_html=True)
    st.markdown('<p class="big-font">Hazard rank is {} because Oxygen Balance is {}.</p>'.format(text, oxy), unsafe_allow_html=True)



def create_dataframe(group, rule_six, oxy, hazard):
    """ To create datafrme and save as csv file

    Args:
        group {[List]} -- A list of all HEFG found in smile
        rule_six {[Integer]} -- value of rule six
        oxy {[Integer]} -- value of oxygen balance
        text {[String]} -- Hazard rank text

    Returns:
        csv file: dataframe saved as csv
    """
    if not group:
        group.append('None')
    df = pd.DataFrame({'HEFG': group, 'Rule Six': rule_six, 'Oxygen Balance': oxy, 'Hazard Rank': hazard})
    return df.to_csv().encode('utf-8')


def main(smile, hefg_list, sym_to_mass, formula):
    """ To create datafrme and save as csv file

    Args:
        smile {[List]} -- A string of Smile.
        hefg_list {[Integer]} -- A list of HEFG to check if present in data.
        sym_to_mass {[Dataframe]} -- Dataframe of atoms mass
        formula {[String]}: string of chemical formula

    Returns:
        hefg {[Integer]} -- number of HEFG present
        r_v {[Integer]} -- value of rule six
        oxy {[Integer]} -- value of oxygen balance
        group {[]} -- A list of all HEFG found in data.
    """
    if smile is not None and formula is not None:
        hefg, group = find_part(smile, hefg_list)
        r_v = rule_six(smile, hefg)
        oxy = oxygen_balance(sym_to_mass, formula)

    return hefg, r_v, oxy, group