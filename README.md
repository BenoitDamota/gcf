# GCF

Generic Cyclic Features (denoted as GCF) are descriptors for molecules focused on cycles. Using the Networkx python library, the vertices (bonds) that do not belong to a cycle are deleted. RDkit is then used to compute the Murcko scaffold on each remaining subgraph that includes a cycle. At this point, the cyclic feature contains the bond and atom types information. In order to work with more generic cyclic features, all atoms with a coordination number of 4 or less are converted to carbon atoms. Since hypervalent carbon produce RDkit errors, hypervalent atoms are left unchanged.

![figure](https://raw.githubusercontent.com/BenoitDamota/gcf/main/gcf_ribavirin.png)

## dependencies

- Numpy
- RDkit 
- Networkx

## usage

```python
aspirin = "O=C(C)Oc1ccccc1C(=O)O"
compute_generic_cyclic_features_with_insat(aspirin)
>>> ['C1=CC=CC=C1']

sitagliptin = "C1CN2C(=NN=C2C(F)(F)F)CN1C(=O)CC(CC3=CC(=C(C=C3F)F)F)N"
compute_generic_cyclic_features_with_insat(sitagliptin)
>>> ['C1=CC=CC=C1', 'C1=CC2CCCCC2=C1']

ribavirin = "C1=NC(=NN1C2C(C(C(O2)CO)O)O)C(=O)N"
compute_generic_cyclic_features_with_insat(ribavirin)
>>> ['C1=CCC=C1', 'C1CCCC1']

```
