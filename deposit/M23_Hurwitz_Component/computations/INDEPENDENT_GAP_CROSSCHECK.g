# Independent character-theoretic cross-check for M23.
# Requires GAP with the Character Table Library (CTblLib).
# The table is the ATLAS-derived library table CharacterTable("M23").
# This script does not import or call any Python/C++ verifier in this release.

LoadPackage("ctbllib");
tbl := CharacterTable("M23");;
names := AtlasClassNames(tbl);;
sizes := SizesConjugacyClasses(tbl);;
irr := Irr(tbl);;

ClassPos := function(name)
  local p;
  p := Position(names,name);
  if p = fail then Error("missing M23 class ",name); fi;
  return p;
end;

# Number of ordered tuples x_i in the listed classes with product equal to a
# fixed element of targetClass.  This is the Frobenius character formula.
FixedProductCount := function(classNames,targetClass)
  local pos,tpos,r,prodSizes,s,chi,val;
  pos := List(classNames,ClassPos);
  tpos := ClassPos(targetClass);
  r := Length(pos);
  prodSizes := Product(List(pos,i->sizes[i]));
  s := 0;
  for chi in irr do
    val := Product(List(pos,i->chi[i]));
    s := s + val*ComplexConjugate(chi[tpos])/(chi[1]^(r-1));
  od;
  return prodSizes*s/Size(tbl);
end;

# Number of ordered product-one tuples with entries in the listed classes.
TupleCount := function(classNames)
  local pos,r,prodSizes,s,chi;
  pos := List(classNames,ClassPos);
  r := Length(pos);
  prodSizes := Product(List(pos,i->sizes[i]));
  s := Sum(irr,chi->Product(List(pos,i->chi[i]))/(chi[1]^(r-2)));
  return prodSizes*s/Size(tbl);
end;

if Size(tbl) <> 10200960 then Error("M23 order mismatch"); fi;
if FixedProductCount(["2A","2A"],"2A") <> 98 then
  Error("2A*2A -> 2A Frobenius count mismatch");
fi;
if FixedProductCount(["3A","3A"],"2A") <> 2688 then
  Error("3A*3A -> 2A Frobenius count mismatch");
fi;
if FixedProductCount(["2A","2A","3A"],"5A") <> 127200 then
  Error("fixed-5A four-branch Frobenius count mismatch");
fi;

Print("M23_GAP_FROBENIUS_CROSSCHECK\n");
Print("PASS True\n");
Print("2A2A_TO_2A 98 = 84+14\n");
Print("3A3A_TO_2A 2688 = 1344+896+224+224\n");
Print("FIXED_5A_2A2A3A 127200\n");
QUIT_GAP(0);
