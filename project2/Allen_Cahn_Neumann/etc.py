import yaml

def read_yaml(file: str):
    # read YAML file
    with open(file, 'r') as f:
        try:
            out = yaml.safe_load(f)
            return out
        except yaml.YAMLError as exc:
            print(exc)

def sub2ind ( array_shape , rows , cols ):
    return rows * array_shape [1] + cols