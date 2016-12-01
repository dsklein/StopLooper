import ROOT

infile_corr = open( "limits_corridor.txt" )
infile_base = open( "limits_baseline.txt" )

h_limits_corr = ROOT.TH2D("limits_corridor", "limits corridor", 37, 87.5, 1012.5, 19, -12.5, 462.5 )
h_limits_base = ROOT.TH2D("limits_baseline", "limits baseline", 37, 87.5, 1012.5, 19, -12.5, 462.5 )


count = 0
for line in infile_corr.readlines():
	fields = line.split()
	limit  = float(fields[4])
	filename = fields[0].split('_')
	m_stop = float(filename[4]) # get the stop mass
	m_lsp  = float(filename[5].split('.')[0]) # get the LSP mass
	binnum = h_limits_corr.FindBin( m_stop, m_lsp )
	h_limits_corr.SetBinContent( binnum, limit )
	count += 1
print 'Found', count, 'mass points for the corridor limits.'

count = 0
for line in infile_base.readlines():
	fields = line.split()
	limit  = float(fields[4])
	filename = fields[0].split('_')
	m_stop = float(filename[4]) # get the stop mass
	m_lsp  = float(filename[5].split('.')[0]) # get the LSP mass
	binnum = h_limits_base.FindBin( m_stop, m_lsp )
	h_limits_base.SetBinContent( binnum, limit )
	count += 1
print 'Found', count, 'mass points for the baseline limits.'

outfile = ROOT.TFile("limit_plots.root", "RECREATE")
outfile.cd()
h_limits_corr.Write()
h_limits_base.Write()
outfile.Close()
