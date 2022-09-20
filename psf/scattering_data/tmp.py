
with open('xray_scattering_params.py','r') as f_in:
    lines = f_in.readlines()

_lines = [x.strip() for x in lines]

nlines = len(lines)
nelem = nlines//2

print(_lines[0])
print(_lines[1])

lines = []

shift = 0
for ii in range(nelem):
    _ = _lines[shift].strip()
    _ += _lines[shift+1].strip()

    _ = [x.strip().strip('}').strip('{') for x in _.split(',')[:-1]]
    lines.append(_)
    
    shift += 2
    
with open('_xrays.py','w') as f_out:
    f_out.write('scattering_params = { \\\n')
    for line in lines:
        elem = line[0].split(':')[1].strip()
        params = line[1:]
        print(params)
        
        _ = f'    {elem:>8} : '+'{' 
        for ii in range(8):
            param = params[ii]
            param = param.split(':')
            label = param[0]
            val = float(param[1])
            _ += f' {label:5} : {val:9.6f},\n                '
        param = params[-1]
        param = param.split(':')
        label = param[0]
        val = float(param[1])
        _ += f' {label:5} : {val:9.6f}'+'},\n'
    
        f_out.write(_)

