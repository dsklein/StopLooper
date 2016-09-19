import ROOT

infile1 = open( "limits_corridor.txt" )
infile2 = open( "limits_ichep.txt" )

h_limits1 = ROOT.TH2D("limits1", "limits1", 37, 87.5, 1012.5, 19, -12.5, 462.5 )
h_limits2 = ROOT.TH2D("limits2", "limits2", 37, 87.5, 1012.5, 19, -12.5, 462.5 )


for line in infile1.readlines():
    fields = line.split()
    limit  = float(fields[4])
    filename = fields[0].split('_')
    m_stop = float(filename[4]) # get the stop mass
    m_lsp  = float(filename[5].split('.')[0]) # get the LSP mass
    binnum = h_limits1.FindBin( m_stop, m_lsp )
    h_limits1.SetBinContent( binnum, limit )


for line in infile2.readlines():
    fields = line.split()
    limit  = float(fields[4])
    filename = fields[0].split('_')
    m_stop = float(filename[4]) # get the stop mass
    m_lsp  = float(filename[5].split('.')[0]) # get the LSP mass
    binnum = h_limits2.FindBin( m_stop, m_lsp )
    h_limits2.SetBinContent( binnum, limit )

outfile = ROOT.TFile("limit_plots.root", "RECREATE")
outfile.cd()
h_limits1.Write()
h_limits2.Write()
outfile.Close()
