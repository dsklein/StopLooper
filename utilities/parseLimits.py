import ROOT

infile_corr  = open( "limits_corridor.txt" )
infile_ichep = open( "limits_ichep.txt" )

h_limits_corr  = ROOT.TH2D("limits_corr",  "limits corridor", 37, 87.5, 1012.5, 19, -12.5, 462.5 )
h_limits_ichep = ROOT.TH2D("limits_ichep", "limits ichep",    37, 87.5, 1012.5, 19, -12.5, 462.5 )


for line in infile_corr.readlines():
    fields = line.split()
    limit  = float(fields[4])
    filename = fields[0].split('_')
    m_stop = float(filename[4]) # get the stop mass
    m_lsp  = float(filename[5].split('.')[0]) # get the LSP mass
    binnum = h_limits_corr.FindBin( m_stop, m_lsp )
    h_limits_corr.SetBinContent( binnum, limit )


for line in infile_ichep.readlines():
    fields = line.split()
    limit  = float(fields[4])
    filename = fields[0].split('_')
    m_stop = float(filename[4]) # get the stop mass
    m_lsp  = float(filename[5].split('.')[0]) # get the LSP mass
    binnum = h_limits_ichep.FindBin( m_stop, m_lsp )
    h_limits_ichep.SetBinContent( binnum, limit )

outfile = ROOT.TFile("limit_plots.root", "RECREATE")
outfile.cd()
h_limits_corr.Write()
h_limits_ichep.Write()
outfile.Close()
