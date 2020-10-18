# Makefile for package xps.
MAIN = xps

include Makefile.arch

#------------------------------------------------------------------------------

# StatUtils
MYHDR       = TMLMath.h TStat.h StatUtils.h
MYOBJ       = TMLMath.$(ObjSuf) TStat.$(ObjSuf) StatUtils.$(ObjSuf)
MYSRC       = TMLMath.$(SrcSuf) TStat.$(SrcSuf) StatUtils.$(SrcSuf)

# XPSBase
MYHDR      += XPSBase.h XPSUtils.h 
MYOBJ      += XPSBase.$(ObjSuf) XPSUtils.$(ObjSuf)
MYSRC      += XPSBase.$(SrcSuf) XPSUtils.$(SrcSuf)

# XPSSchemes
MYHDR      += XPSSchemes.h
MYOBJ      += XPSSchemes.$(ObjSuf)
MYSRC      += XPSSchemes.$(SrcSuf)

# XPSDataTypes
MYHDR      += XPSDataTypes.h
MYOBJ      += XPSDataTypes.$(ObjSuf)
MYSRC      += XPSDataTypes.$(SrcSuf)

# XProjectHandler
MYHDR      += XPSProjectHandler.h
MYOBJ      += XPSProjectHandler.$(ObjSuf)
MYSRC      += XPSProjectHandler.$(SrcSuf)

# XPSData
MYHDR      += XPSData.h
MYOBJ      += XPSData.$(ObjSuf)
MYSRC      += XPSData.$(SrcSuf)

# XPSProcessing
MYHDR      += XPSProcessing.h XPSHybridizer.h XPSSelector.h XPSNormalizer.h XPSPreProcessing.h XPSNormation.h
MYOBJ      += XPSProcessing.$(ObjSuf) XPSHybridizer.$(ObjSuf) XPSSelector.$(ObjSuf) XPSNormalizer.$(ObjSuf) XPSPreProcessing.$(ObjSuf) XPSNormation.$(ObjSuf)
MYSRC      += XPSProcessing.$(SrcSuf) XPSHybridizer.$(SrcSuf) XPSSelector.$(SrcSuf) XPSNormalizer.$(SrcSuf) XPSPreProcessing.$(SrcSuf) XPSNormation.$(SrcSuf)

# XPSAnalysis
MYHDR      += XPSFilter.h XPSAnalyzer.h XPSAnalysis.h 
MYOBJ      += XPSFilter.$(ObjSuf) XPSAnalyzer.$(ObjSuf) XPSAnalysis.$(ObjSuf)
MYSRC      += XPSFilter.$(SrcSuf) XPSAnalyzer.$(SrcSuf)XPSAnalysis.$(SrcSuf)

# xps
MYHDR      += rwrapper.h $(MAIN)LinkDef.h
MYOBJ      += rwrapper.$(ObjSuf) $(MAIN)Dict.$(ObjSuf)
MYSRC      += rwrapper.$(SrcSuf) $(MAIN)Dict.$(SrcSuf)

ifeq ($(PLATFORM),win32)
MYLIBS = '$(ROOTSYS)/lib/libRIO.lib $(ROOTSYS)/lib/libGed.lib $(ROOTSYS)/lib/libTreePlayer.lib $(ROOTSYS)/lib/libTreeViewer.lib'  
else
MYLIBS = -lRIO -lGed -lTreePlayer -lTreeViewer
endif

MYDLL       = $(MAIN).$(DllSuf)
ifeq ($(PLATFORM),win32)
MYLIB       = $(MAIN).lib
else
MYLIB       = $(MYDLL)
endif

OBJS        = $(MYOBJ)
PROGRAMS    = $(MYDLL)

#------------------------------------------------------------------------------

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

all:            $(PROGRAMS)

$(MYDLL):    $(MYOBJ)
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
		$(LD) $(SOFLAGS)$@ $(LDFLAGS) $^ $(GLIBS) $(MYLIBS) $(OutPutOpt) $@
ifneq ($(subst $(MACOSX_MINOR),,1234),1234)
ifeq ($(MACOSX_MINOR),4)
		ln -sf $@ $(subst .$(DllSuf),.so,$@)
else
		$(LD) -bundle $(LDFLAGS) $^ $(GLIBS) $(MYLIBS) \
		   $(OutPutOpt) $(subst .$(DllSuf),.so,$@)
endif
endif
else
ifeq ($(PLATFORM),win32)
		bindexplib $* $^ > $*.def
		lib -nologo -MACHINE:IX86 $^ -def:$*.def \
		   $(OutPutOpt)$(MYLIB)
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $*.exp $(LIBS) \
		   $(OutPutOpt)$@
else
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(GLIBS) $(MYLIBS) $(OutPutOpt) $@ 
endif
endif

clean-$(MAIN):
		@rm -f $(OBJS) core $(MAIN)Dict.* $(MAIN).$(DllSuf) *.so

mostlyclean:    clean
clean::         clean-$(MAIN)

.SUFFIXES: .$(SrcSuf)

###


$(MAIN)Dict.$(SrcSuf): $(MYHDR)
	@echo "Generating dictionary $@..."
	@rootcint -f $@ -c $^

.$(SrcSuf).$(ObjSuf):
	$(CXX) -I${ROOTSYS}/include $(CXXFLAGS) -c $<
