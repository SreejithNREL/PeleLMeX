namespace amrex {

const char* buildInfoGetBuildDate() {

  static const char BUILD_DATE[] = __DATE__ " "  __TIME__;
  return BUILD_DATE;
}

const char* buildInfoGetBuildDir() {

  static const char BUILD_DIR[] = "E:/PeleLMeX_Dev/buildHotBubble/Exec/RegTests/HotBubble";
  return BUILD_DIR;
}

const char* buildInfoGetBuildMachine() {

  static const char BUILD_MACHINE[] = "Windows SreejithNA  Professional  (Build 26100) AMD64";
  return BUILD_MACHINE;
}

const char* buildInfoGetAMReXDir() {

  static const char AMREX_DIR[] = "E:/PeleLMeX_Dev/Submodules/PelePhysics/Submodules/amrex";
  return AMREX_DIR;
}

const char* buildInfoGetComp() {

  static const char COMP[] = "MSVC";
  return COMP;
}

const char* buildInfoGetCompVersion() {

  static const char COMP_VERSION[] = "19.42.34435.0";
  return COMP_VERSION;
}

// deprecated
const char* buildInfoGetFcomp() {

  static const char FCOMP[] = "";
  return FCOMP;
}

// deprecated
const char* buildInfoGetFcompVersion() {

  static const char FCOMP_VERSION[] = "";
  return FCOMP_VERSION;
}

const char* buildInfoGetCXXName() {

  static const char CXX_comp_name[] = "C:/Program Files/Microsoft Visual Studio/2022/Community/VC/Tools/MSVC/14.42.34433/bin/Hostx64/x64/cl.exe";
  return CXX_comp_name;
}

const char* buildInfoGetFName() {

  static const char F_comp_name[] = "";
  return F_comp_name;
}

const char* buildInfoGetCXXFlags() {

  static const char CXX_flags[] = "-DUSE_CONSTANT_TRANSPORT -DUSE_FUEGO_EOS -IE:/PeleLMeX_Dev/Submodules/PelePhysics/Source/Utility/TurbInflow -IE:/PeleLMeX_Dev/Submodules/PelePhysics/Source/Utility/TurbForcing -IE:/PeleLMeX_Dev/Submodules/PelePhysics/Source/Utility/Diagnostics -IE:/PeleLMeX_Dev/Submodules/PelePhysics/Source/Utility/PltFileManager -IE:/PeleLMeX_Dev/Submodules/PelePhysics/Source/Utility/BlackBoxFunction -IE:/PeleLMeX_Dev/Submodules/PelePhysics/Source/Utility/PMF -IE:/PeleLMeX_Dev/Submodules/PelePhysics/Source/Utility/Utilities -IE:/PeleLMeX_Dev/Submodules/PelePhysics/Submodules/amrex/Src/Extern/SUNDIALS -IE:/PeleLMeX_Dev/Submodules/PelePhysics/Source -IE:/PeleLMeX_Dev/Submodules/PelePhysics/Source/Transport -IE:/PeleLMeX_Dev/Submodules/PelePhysics/Source/Eos -IE:/PeleLMeX_Dev/Submodules/PelePhysics/Mechanisms/air -IE:/PeleLMeX_Dev/Submodules/PelePhysics/Source/Reactions /O2 /Ob2 /DNDEBUG /DWIN32 /D_WINDOWS /EHsc -W4 -D_XOPEN_SOURCE";
  return CXX_flags;
}

const char* buildInfoGetFFlags() {

  static const char F_flags[] = "";
  return F_flags;
}

const char* buildInfoGetLinkFlags() {

  static const char link_flags[] = "";
  return link_flags;
}

const char* buildInfoGetLibraries() {

  static const char libraries[] = "";
  return libraries;
}

const char* buildInfoGetAux(int i) {

  //static const char AUX1[] = "${AUX[1]}";
  
  static const char EMPT[] = "";

  switch(i)
  {
    
    default: return EMPT;
  }
}

int buildInfoGetNumModules() {
  // int const num_modules = X;
  int const num_modules = 0;
  return num_modules;
}

const char* buildInfoGetModuleName(int i) {

  //static const char MNAME1[] = "${MNAME[1]}";
  
  static const char EMPT[] = "";

  switch(i)
  {
    
    default: return EMPT;
  }
}

const char* buildInfoGetModuleVal(int i) {

  //static const char MVAL1[] = "${MVAL[1]}";
  
  static const char EMPT[] = "";

  switch(i)
  {
    
    default: return EMPT;
  }
}

const char* buildInfoGetGitHash(int i) {

  //static const char HASH1[] = "${GIT[1]}";
  static const char HASH1[] = "v23.05-214-gba5babba4622-dirty";
  static const char HASH2[] = "25.07";
  static const char EMPT[] = "";

  switch(i)
  {
    case 1: return HASH1;
    case 2: return HASH2;
    default: return EMPT;
  }
}

const char* buildInfoGetBuildGitHash() {

  //static const char HASH[] = "${GIT}";
  static const char HASH[] = "";

  return HASH;
}

const char* buildInfoGetBuildGitName() {

  //static const char NAME[] = "";
  static const char NAME[] = "";

  return NAME;
}

}
