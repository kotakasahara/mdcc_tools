#include "FeatureDefinition.h"

FeatureDefinition::FeatureDefinition(int inType, int inColumn){
  type=inType;
  column=inColumn;
}
FeatureDefinition::FeatureDefinition(int inType, int inColumn, string inH){
  type=inType;
  column=inColumn;
  header=inH;
}
