import svg_stack as ss
import os

doc = ss.Document()

layout1 = ss.HBoxLayout()
layout1.addSVG('fig_modelPredictionBurstsIndividual_VenancesBin0.svg',alignment=ss.AlignTop|ss.AlignHCenter)
#layout1.addSVGNoLayout('protocolIllustration.svg',x=0,y=-100)
layout1.addSVGNoLayout('spike-trains_illustration.svg',x=-510,y=-105)
#layout1.addSVGNoLayout('fig_modelPredictionBurstsIndividual_VenancesBin0.svg',x=00,y=00)
#layout1.addSVGNoLayout('spatial_illustration.svg',x=-1450,y=-35)

#layout2 = ss.VBoxLayout()

#layout2.addSVG('red_ball.svg',alignment=ss.AlignCenter)
#layout2.addSVG('red_ball.svg',alignment=ss.AlignCenter)
#layout2.addSVG('red_ball.svg',alignment=ss.AlignCenter)
#layout1.addLayout(layout2)

doc.setLayout(layout1)

figname = 'fig_modelPredictionBurstsIndividual_VenancesBin0_complete'

doc.save(figname+'.svg')
os.system('inkscape -f '+str(figname)+'.svg -A '+str(figname)+'.pdf')
os.system('convert '+str(figname)+'.pdf -quality 100 '+str(figname)+'.png')
