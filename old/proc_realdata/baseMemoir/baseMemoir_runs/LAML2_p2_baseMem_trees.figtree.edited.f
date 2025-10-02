#NEXUS
begin taxa;
	dimensions ntax=30;
	taxlabels
	1
	10
	11
	12
	13[&!color=#ff0000]
	14
	15[&!color=#ff0000]
	16
	17[&!color=#ff0000]
	18[&!color=#ff0000]
	20
	23[&!color=#ff0000]
	24
	25[&!color=#ff0000]
	26[&!color=#ff0000]
	27[&!color=#ff0000]
	28[&!color=#ff0000]
	29
	3
	30
	31
	32
	33
	34
	35
	4
	5
	6
	8[&!color=#ff0000]
	9
;
end;

begin trees;
	tree tree_1 = [&R] ((((((((12:0.063482,3:0.063482):0.441719,(13:0.006008,6:0.006008):0.499193):0.421133,((25:0.669226,16:0.669226):0.232623,((24:0.485662,30:0.485662):0.412538,(((31:0.01006,4:0.01006):0.472564,1:0.482624):0.274856,23:0.75748):0.14072):0.003649):0.024485):5.95E-4,(((17:0.143886,((26:4.7E-4,28:4.7E-4):4.7E-4,27:9.39E-4):0.142947):0.425966,8:0.569852):0.325877,29:0.895729):0.031199):0.011405,34:0.938333):0.021649,((20:0.662931,11:0.662931):0.225834,(((18:0.427357,32:0.427357):0.411225,(33:0.194535,35:0.194535):0.644047):0.04829,(5:0.202902,10:0.202902):0.683971):0.001893):0.071216):0.017515,((14:0.232887,9:0.232887):0.19554,15:0.428426):0.549071):0.022503);
end;

begin figtree;
	set appearance.backgroundColorAttribute="Default";
	set appearance.backgroundColour=#ffffff;
	set appearance.branchColorAttribute="User selection";
	set appearance.branchColorGradient=false;
	set appearance.branchLineWidth=1.0;
	set appearance.branchMinLineWidth=0.0;
	set appearance.branchWidthAttribute="Fixed";
	set appearance.foregroundColour=#000000;
	set appearance.hilightingGradient=false;
	set appearance.selectionColour=#2d3680;
	set branchLabels.colorAttribute="User selection";
	set branchLabels.displayAttribute="Branch times";
	set branchLabels.fontName="sansserif";
	set branchLabels.fontSize=8;
	set branchLabels.fontStyle=0;
	set branchLabels.isShown=false;
	set branchLabels.significantDigits=4;
	set layout.expansion=0;
	set layout.layoutType="RADIAL";
	set layout.zoom=0;
	set legend.attribute=null;
	set legend.fontSize=10.0;
	set legend.isShown=false;
	set legend.significantDigits=4;
	set nodeBars.barWidth=4.0;
	set nodeBars.displayAttribute=null;
	set nodeBars.isShown=false;
	set nodeLabels.colorAttribute="User selection";
	set nodeLabels.displayAttribute="Node ages";
	set nodeLabels.fontName="sansserif";
	set nodeLabels.fontSize=8;
	set nodeLabels.fontStyle=0;
	set nodeLabels.isShown=false;
	set nodeLabels.significantDigits=4;
	set nodeShapeExternal.colourAttribute="User selection";
	set nodeShapeExternal.isShown=true;
	set nodeShapeExternal.minSize=10.0;
	set nodeShapeExternal.scaleType=Width;
	set nodeShapeExternal.shapeType=Circle;
	set nodeShapeExternal.size=4.0;
	set nodeShapeExternal.sizeAttribute="Fixed";
	set nodeShapeInternal.colourAttribute="User selection";
	set nodeShapeInternal.isShown=false;
	set nodeShapeInternal.minSize=10.0;
	set nodeShapeInternal.scaleType=Width;
	set nodeShapeInternal.shapeType=Circle;
	set nodeShapeInternal.size=4.0;
	set nodeShapeInternal.sizeAttribute="Fixed";
	set polarLayout.alignTipLabels=false;
	set polarLayout.angularRange=0;
	set polarLayout.rootAngle=0;
	set polarLayout.rootLength=100;
	set polarLayout.showRoot=true;
	set radialLayout.spread=0.0;
	set rectilinearLayout.alignTipLabels=false;
	set rectilinearLayout.curvature=0;
	set rectilinearLayout.rootLength=100;
	set scale.offsetAge=0.0;
	set scale.rootAge=1.0;
	set scale.scaleFactor=1.0;
	set scale.scaleRoot=false;
	set scaleAxis.automaticScale=true;
	set scaleAxis.fontSize=8.0;
	set scaleAxis.isShown=false;
	set scaleAxis.lineWidth=1.0;
	set scaleAxis.majorTicks=0.25;
	set scaleAxis.minorTicks=0.05;
	set scaleAxis.origin=0.0;
	set scaleAxis.reverseAxis=false;
	set scaleAxis.showGrid=true;
	set scaleBar.automaticScale=true;
	set scaleBar.fontSize=10.0;
	set scaleBar.isShown=false;
	set scaleBar.lineWidth=1.0;
	set scaleBar.scaleRange=0.2;
	set tipLabels.colorAttribute="User selection";
	set tipLabels.displayAttribute="Names";
	set tipLabels.fontName="sansserif";
	set tipLabels.fontSize=20;
	set tipLabels.fontStyle=0;
	set tipLabels.isShown=true;
	set tipLabels.significantDigits=4;
	set trees.order=false;
	set trees.orderType="increasing";
	set trees.rooting=false;
	set trees.rootingType="User Selection";
	set trees.transform=false;
	set trees.transformType="cladogram";
end;

