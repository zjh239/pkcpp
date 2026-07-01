This is a C++ migration of the Polikit atomistic analysis package.

## Migration progress

✅ RDF
✅ CN
❌ Poly-analysis
❌ BAD
❌ RSA 

## Plan

I realize that I forgot what I did with the code after several months. So I plan to give the code ability to read TOML file as analysis input script.

## Howto

### Build

Starting from directory of the code, there are two ways to build: 

- **cmake (3.29+) & Ninja (1.13.0+)**

```bash
mkdir build && cd build
cmake -GNinja ../
```

- **xmake (3.0.0+)**

```bash
xmake
```

### Usage

```bash
polikit -f [data] -p [pbc] -*
```

`-*` means computing options, they include

```
-nf [cutoffs] # neighbor finding.
-rdf [cutoffs]
-cn [cutoffs]
-d2min [cutoffs]
-poly [cutoffs]
-bad [cutoffs]
-ring [cutoffs]
-cluster [cutoffs]
```
#### Example

`polikit -f ../test/ga2o3.xyz -p 1`
