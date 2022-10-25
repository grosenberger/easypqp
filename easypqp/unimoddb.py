import click
from tqdm import tqdm

# Unimod parsing
import xml.etree.cElementTree as ET


def site_validation(site_input):
    """
    Perform a check to ensure inputs are valid
    Arguments:
        site_input: (list) list of amino acid residues, or terminal notation, or wild card notation (*).
    Returns:
        Nothing is returned. An error is raised if the input contains a non-valid site.
    """
    acceptable_sites = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", 'U', 'O', '[', ']', 'n', 'c', '*']
    site_check = [site not in acceptable_sites for site in site_input]
    if any(site_check):
        raise click.ClickException( f"Incorrect site specificity input, site(s) {', '.join([i for (i, v) in zip(site_input, site_check) if v])} is not valid. Acceptable sites: {', '.join(acceptable_sites)}")

def site_specificity_transform(site_input):
    """
    Transform input site to return the site and position. Transforms terminal notation to site notation in unimod.xml and whether its any terminal site or a protein terminal site.
    Arguments:
        site_input: (list) list of amino acid residues, or terminal notation, or wild card notation (*).
    Returns:
        Returns a tuple of list of sites and list of positions
    """
    # Site and Position Mapping
    terminal_map = {'[':'N-term', ']':'C-term', 'n':'N-term', 'c':'C-term'}
    site_position_map = {'[':'Protein N-term', ']':'Protein C-term', 'n':'Any N-term', 'c':'Any C-term'}
    # Split sites
    site_input = [site for site in site_input]
    site_validation(site_input)
    sites=[]; positions=[]
    for site in site_input:
        if site in terminal_map.keys():
            sites.append(terminal_map[site])
            positions.append(site_position_map[site])
        elif site=="*":
            sites.append("*")
            positions.append("*")
        else:
            sites.append(site)
            positions.append("Anywhere")
    return sites, positions

def unimod_filter(unimod_file, out_file, accession_ids, site_specificity):
    """
    Filter an input unimod to restrict for specific modifications and site specificities
    Arguments:
        unimod_file: (str) path/filename of input unimod.xml file.
        out_file: (str) path/filename to write out new filtered unimod.xml file
        accession_ids: (list) list of unimod accession ids to restrict for. i.e. ['1','21','35]
        site_specificity: (list) list of site specificties to further restrict corresponding unimod for. i.e. ['n','STY','M], will restrict acetylation for any N-Term, phosphorylation for serine, threonine, and tyrosine, and oxidation for methionine.
    Returns:
        Nothing is returned. The restricted unimod database is written to the out_file.
    """
    # Register Namespace
    ET.register_namespace('umod', 'http://www.unimod.org/xmlns/schema/unimod_2')

    # Read in unimod XML database
    click.echo(f"INFO: Loading XML data from {unimod_file}")
    tree = ET.parse(unimod_file)
    root = tree.getroot()

    # Namespace
    ns = {'umod':'http://www.unimod.org/xmlns/schema/unimod_2'}

    # Generate root for new filtered unimod XML
    root_out = ET.Element(root.tag, root.attrib)

    # Append elements subelements
    umod_elements = root.findall("umod:elements", ns)
    root_out.append(umod_elements[0])

    # Append desired modifications
    mod_entries = root.findall('umod:modifications', ns)[0]
    mod_out = ET.Element(mod_entries.tag, mod_entries.attrib)
    i=0
    pbar = tqdm(accession_ids)
    pbar_desc = "INFO: Restricting"
    for record_id in pbar:
        add_unimod_entry = mod_entries.findall(f"./umod:mod/[@record_id='{record_id}']", ns)[0]
        if site_specificity is not None:
            site, position = site_specificity_transform(site_specificity[i])
            # Update progess bar description
            pbar_desc = f"INFO: Restricting..{add_unimod_entry.attrib.get('title')}({','.join(site)})"
            pbar.set_description(pbar_desc)
            if site != "*":
                for unimod_site in add_unimod_entry.findall(f"./umod:specificity", ns):
                    if unimod_site.attrib['site'] in site and unimod_site.attrib['position'] in position:
                        # If current specificity element is a requested one, continue on
                        continue
                    else:
                        # Remove specificities that do not match requested specificities
                        add_unimod_entry.remove(unimod_site)
        else:
            # Update progess bar description
            pbar_desc = f"INFO: Restricting..{add_unimod_entry.attrib.get('title')}"
            pbar.set_description(pbar_desc)
        # click.echo(f"INFO: Appending to filtered unimod XML - title={add_unimod_entry.attrib.get('title')} with record_id={add_unimod_entry.attrib.get('record_id')}")
        mod_out.append( add_unimod_entry )
        i+=1
    root_out.append(mod_out)

    # Append amino acids
    umod_amino_acids = root.findall("umod:amino_acids", ns)
    root_out.append(umod_amino_acids[0])

    # Append mod bricks
    umod_mod_bricks = root.findall("umod:mod_bricks", ns)
    root_out.append(umod_mod_bricks[0])

    # Generate element hierarchy to write out to xml
    tree_out = ET.ElementTree(root_out)
    # For Pretty-Printing
    ET.indent(tree_out, '  ')
    # Write out filtered unimod xml database
    click.echo(f"INFO: Writing out filtered unimod XML file to {out_file}")
    tree_out.write(out_file, encoding="UTF-8", xml_declaration=True, method="xml")

    # Insert Top Comment
    # TODO: This may not be the best way to add the top level comment in standard unimod.xml database files. Might be able to use lxml instead, requiring an additional dependency
    with open(out_file, 'r+', encoding="utf-8") as file_handle:
        lines = file_handle.readlines()     
        lines.insert(1, "<!--Copyright (C) 2002-2006 Unimod; this information may be copied, distributed and/or-->\n<!--modified under certain conditions, but it comes WITHOUT ANY WARRANTY; see the-->\n<!--accompanying Design Science License for more details-->\n")  # you can use any index if you know the line index
        file_handle.seek(0)                 
        file_handle.writelines(lines)       