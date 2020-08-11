from augur.utils import read_metadata
import sys
import json

fields_to_add = ['coverage', 'date_seq', 'lab']

if __name__ == "__main__":
    data = {'nodes': {}}
    ms_dict, ms_columns = read_metadata(sys.argv[1])
    private_dict, _ = read_metadata(sys.argv[2])

    print("\t".join(ms_columns+fields_to_add))
    
    for strain, data in ms_dict.items():
        line = [data[f] for f in ms_columns]
        for key in fields_to_add:
            try:
                line.append(str(private_dict[strain][key]))
            except KeyError:
                line.append("")
        print("\t".join(line))

