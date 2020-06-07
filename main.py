import argparse
import json
import math


def parse_init():
    parser = argparse.ArgumentParser(description='Discrete source model.')
    parser.add_argument(
        '-f', '--file',
        type=argparse.FileType('r', encoding='UTF-8'),
        required=True,
        help='Name of the discrete source description file. Required.'
    )
    parser.add_argument(
        '-r',
        type=float,
        nargs='?',
        default=None,
        help='Coding speed R.'
    )
    parser.add_argument(
        '-e', '--epsilon',
        type=float,
        nargs='?',
        help='Decoding error bound Epsilon.'
    )
    parser.add_argument(
        '-q',
        type=int,
        nargs='?',
        help="Power q decoder's alphabet."
    )
    args = parser.parse_args()
    #print(args)
    if args.file:
        if args.file.name.split('.')[-1]!='json':
            print("Warning: your file doesn't have 'json' extension. Continue?[yes/no]")
            input_flag = True
            while input_flag:
                answer = input().lower()
                if answer in ['n', 'no', 'нет', 'н']:
                    args.file.close()
                    exit(0)
                elif answer in ['y', 'yes', 'да', 'д']:
                    input_flag = False
                else:
                    print("Please type [yes/no]")
    #print(args.r, args.q, args.epsilon)
    if args.r and args.q and args.epsilon:
        return args, True
    return args, False


def check_distribution(data):
    for key, value in data.items():
        prob = 0
        for probability in data[key].values():
            prob += probability
            if probability > 1.0 or probability < 0.0:
                return False
        if prob!=1.0:
            return False
    return True


def check_names_dict(names, data):
    for values in data.values():
        for namekey in values.keys():
            if namekey not in names:
                return False
    return True


def check_names_list(names, list):
    for elem in list:
        if elem not in names:
            return False
    return True


def string_float(s):
    if '/' in s:
        return string_float(s.split('/')[0])/string_float(s.split('/')[1])
    else:
        return float(s)


def get_json(file):
    data = json.load(file)
    file.close()
    keys = data.keys()
    if 'models' not in keys or\
        'switches' not in keys or\
        'source' not in keys:
        exit("Error: json file requires 'models', 'switches' and 'source' keys")
    models = data['models']
    switches = data['switches']
    for key in models.keys():
        for key_m, value_m in models[key].items():
            models[key][key_m] = string_float(value_m)
    for key in switches.keys():
        for key_m, value_m in switches[key].items():
            switches[key][key_m] = string_float(value_m)
    if not check_distribution(data['models']) or not check_distribution(data['switches']):
        exit("Error: wrong source distribution function.")
    if not check_names_dict(models.keys(), data['switches']) or not check_names_list(switches.keys(), [i['switch'] for i in data['source']]):
        exit("Error: unknown name is used.")
    return data


def first_mode(data):
    for source in data['source']:
        if source['input'] != []:
            exit("Error: non stationary zero-memory source is provided for first mode.")
    switch = data['source'][0]['switch']
    alphabet = dict()
    size = 0
    for key, value in data['switches'][switch].items():
        probs = data['models'][key]
        for symb, prob_value in probs.items():
            if symb not in alphabet:
                alphabet[symb] = {"num": size, "p": value * prob_value}
                size += 1
            else:
                alphabet[symb]['p'] += value * prob_value
            if alphabet[symb]['p'] > 1.0:
                exit("Error: Wrong probability in model.")
    entropy = -sum(value['p'] * math.log2(value['p']) for value in alphabet.values())
    return entropy, alphabet, size


def information(alphabet):
    info = dict()
    for symb, prob in alphabet.items():
        info[symb] = -math.log2(prob['p'])
    return info


def dispersion(alphabet, h):
    mi2 = 0
    for key, value in alphabet.items():
        mi2 += value['p'] * (math.log2(value['p']) ** 2)
    return mi2 - h ** 2

# def get_max_memory(source):
#    for elem in source:
#        for item in elem['input']:


def info_block(l, i, n):
    info = 0
    for t in l:
        info += i[t[0]]
    return info / n



def add_list(l, q, a, n):
    i = n - 1
    while i >= 0:
        if l[i][1] < q - 1:
            l[i] = (next(key for key, value in a.items() if value['num'] == l[i][1] + 1), l[i][1] + 1)
            break
        else:
            l[i] = (next(key for key, value in a.items() if value['num'] == 0), 0)
        i -= 1
    return l


def set_prob(coded_set, alphabet):
    res = 0
    for item in coded_set:
        prob = 1
        for symb in item:
            prob = prob * alphabet[symb]['p']
        res += prob
    return res


def encode(n, delta, h, size, alphabet, info):
    zero_elem = next(key for key, value in alphabet.items() if value['num'] == 0)
    list_j = [(zero_elem, 0)] * n
    coded_set = list()
    for j in range(size ** n):
        # print("done")
        list_j = add_list(list_j, size, alphabet, n)
        j_info = info_block(list_j, info, n)
        if h - delta <= j_info <= h + delta:
            coded_set.append([elem[0] for elem in list_j])
    return coded_set


def second_mode(args, data):
    h, alphabet, size = first_mode(data)
    print("source entropy:", h)
    if args.r < h:
        exit("Error: cannot build required coder: R < H")
    size = len(alphabet)
    print("alphabet size:", size)
    info = information(alphabet)
    d = dispersion(alphabet, h)
    print("Dispersion:", d)
    temp = 1 / ((args.r - h) ** 2) * (d / args.epsilon)
    if int(temp) == temp:
        n_min = int(temp)
    else:
        n_min = int(temp) + 1
    n = None
    if n_min > int(args.q / h):
        print("Cannot find necessary value in base calculation. Perform iteration.")
        temp = args.q / args.r
        if int(temp) == temp:
            low = int(temp)
        else:
            low = int(temp) + 1
        for i in range(low, int(args.q / h) + 1):
            delta = (d / (i * args.epsilon)) ** 0.5
            coded_set = encode(i, delta, h, size, alphabet, info)
            true_eps = 1 - set_prob(coded_set, alphabet)
            if len(coded_set) <= 2 ** args.q and true_eps <= args.epsilon:
                n = i
                break
            if len(coded_set) > 2 ** args.q and true_eps <= args.epsilon:
                print("n = {}; unable to achieve condition for code. Selecting smaller delta.".format(i))
                #print(delta)
                delta = delta / 2
                #print(delta)
                coded_set = encode(i, delta, h, size, alphabet, info)
                true_eps = 1 - set_prob(coded_set, alphabet)
                #print(len(coded_set))
                if len(coded_set) <= 2 ** args.q and true_eps <= args.epsilon:
                    n = i
                    break
            elif true_eps > args.epsilon and len(coded_set) <= 2 ** args.q:
                delta = delta + (1 - delta) / size
                coded_set = encode(i, delta, h, size, alphabet, info)
                true_eps = 1 - set_prob(coded_set, alphabet)
                print(len(coded_set))
                if len(coded_set) <= 2 ** args.q and true_eps <= args.epsilon:
                    n = i
                    break
        if n is None:
            exit("Error: Cannot build the coder. No available n for required alphabet power.")
        print("Found n:", n)
        print("true epsilon:", true_eps)
        print("coded set length:", len(coded_set))
        form = "{" + "0:0{}b".format(args.q) + "}"
        i = 0
        with open(args.file.name.rsplit('.', 1)[0] + '.code', 'w', encoding='UTF-8') as f:
            for block in coded_set:
                f.write('[' + " ".join(block) + ']    ' + form.format(i) + '\n')
                i += 1
        return 0
    temp = args.q / args.r
    if int(temp) == temp:
        lower_bound = int(temp)
    else:
        lower_bound = int(temp) + 1
    for i in range(max(n_min, lower_bound), int(args.q / h) + 1):
        print(i)
        delta = (d / (i * args.epsilon)) ** 0.5
        coded_set = encode(i, delta, h, size, alphabet, info)
        true_eps = 1 - set_prob(coded_set, alphabet)
        if len(coded_set) <= 2 ** args.q:
            n = i
            break
    if n is None:
        exit("Error: Cannot build coder with required alphabet power")
    print("Found n:", n)
    print("coded set length:", len(coded_set))
    print("true epsilon:", true_eps)
    form = "{" + "0:0{}b".format(args.q) + "}"
    i = 0
    with open(args.file.name.rsplit('.', 1)[0] + '.code', 'w', encoding='UTF-8') as f:
        for block in coded_set:
            f.write('[' + " ".join(block) + ']    ' + form.format(i) + '\n')
            i += 1
    return 0


def test(args, data):
    h, alphabet, size = first_mode(data)
    if args.r < h:
        exit("Error: cannot build required coder: R < H")
    size = len(alphabet)
    info = information(alphabet)
    d = dispersion(alphabet, h)
    print(args.r, h)
    temp = 1 / ((args.r - h) ** 2) * (d / args.epsilon)
    if int(temp) == temp:
        n_min = int(temp)
    else:
        n_min = int(temp) + 1
    print(n_min)
    n = None
    for i in range(12, 21):
        delta = (d / (i * args.epsilon)) ** 0.5
        coded_set = encode(i, delta, h, size, alphabet, info)
        if len(coded_set) <= 2 ** args.q:
            n = i
            break
    if n is None:
        exit("Error: Cannot build coder with required alphabet power")
    print("Found n:", n)
    print("coded set:", len(coded_set))
    print("true epsilon:", 1 - set_prob(coded_set, alphabet))
    form = "{" + "0:0{}b".format(args.q) + "}"
    i = 0
    with open('./codes/' + args.file.name.rsplit('.', 1)[0] + '.code', 'w', encoding='UTF-8') as f:
        for block in coded_set:
            f.write('[' + " ".join(block) + ']    ' + form.format(i) + '\n')
            i += 1
    return 0

def main():
    args, mode = parse_init()
    print("Processing file...")
    data = get_json(args.file)
    print("Done.")
    if mode:
        print("Executing second mode..")
        #test(args, data)
        second_mode(args, data)
    else:
        print("Executing first mode..")
        print("Stationary zero-memory source entropy is:", first_mode(data)[0])
    return 0



if __name__ == "__main__":
    main()

'''form = "{" + "0:0{}b".format(15) + "}"
i = 143
print(form.format(i))'''