# Given a private DRC EBOV auspice JSON, convert the JSON to one with
# names matching the public build
# Name format: Private name: lab-<LABID>[_epi-<EPIID>]. public name: <LABID>
#
# Argument 1: private JSON path. v2 auspice format _or_ v1 tree JSON
# prints JSON to STDOUT. Same format as input


import json
import sys
import re

name_regex = r'lab-([-A-Za-z0-9]+)(_epi-([-A-Za-z0-9]+))?$'

## recursive function to change name
def change_node_name(node):
  v1 = 'strain' in node.keys()
  m = re.match(name_regex, node['strain' if v1 else 'name'])
  if m:
    # print(f"{node['name']} -> {m.groups()[0]}")
    node['strain' if v1 else 'name'] = m.groups()[0] ## lab ID
  # else:
        # print(f"{node['name']} unchanged")

  if 'children' in node.keys():
    for child in node['children']:
      change_node_name(child)

if __name__ == "__main__":
  with open(sys.argv[1]) as fh:
    private = json.load(fh)
  change_node_name(private["tree"])
  print(json.dumps(private, indent=2))