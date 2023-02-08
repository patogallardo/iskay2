#! /usr/bin/env python
from iskay2 import paramfile
import json

p = paramfile.load_paramfile()

queries = ["lum > 4.3e10",
           "lum > 6.1e10",
           "lum > 7.9e10",
           "lum > 4.3e10 and lum < 6.1e10",
           "lum > 6.1e10 and lum < 7.9e10",
           "lum > 4.3e10 and z<0.5",
           "lum > 4.3e10 and z>0.5",
           "lum > 4.3e10 and z>0.44 and z<0.66"]

base_query = p["QUERY"]
base_name = p["NAME"]
for query in queries:
    newquery = base_query + " and " + query
    new_name = base_name + "_" + query.replace(" ", '_').replace(">", "gt").replace("<", 'lt').replace(".", 'p')
    p["QUERY"] = newquery
    p["NAME"] = new_name
    fname_out = "params_" + query.replace(" ", '_').replace(">", "gt").replace("<", 'lt').replace(".", 'p') + ".json"
    with open(fname_out, 'w') as f:
        json.dump(p, f, indent=4)
