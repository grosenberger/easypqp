import click

# Unimod parsing
import xml.etree.cElementTree as ET
from xml.etree.cElementTree import iterparse

def unimod_filter(unimod_file, out_file, accession_ids):
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
    for record_id in accession_ids:
        # print(record_id)
        add_unimod_entry = mod_entries.findall(f"./umod:mod/[@record_id='{record_id}']", ns)[0]
        click.echo(f"INFO: Appending to filtered unimod XML - title={add_unimod_entry.attrib.get('title')} with record_id={add_unimod_entry.attrib.get('record_id')}")
        mod_out.append( add_unimod_entry )
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
    with open(out_file, 'r+') as file_handle:
        lines = file_handle.readlines()     
        lines.insert(1, "<!--Copyright (C) 2002-2006 Unimod; this information may be copied, distributed and/or-->\n<!--modified under certain conditions, but it comes WITHOUT ANY WARRANTY; see the-->\n<!--accompanying Design Science License for more details-->\n")  # you can use any index if you know the line index
        file_handle.seek(0)                 
        file_handle.writelines(lines)       