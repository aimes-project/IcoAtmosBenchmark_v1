// This file is part of MD-REAL-IO.

//
// MD-REAL-IO is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MD-REAL-IO is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with MD-REAL-IO.  If not, see <http://www.gnu.org/licenses/>.
//
// Author: Julian Kunkel
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "option.h"
static int print_value(option_help * o)
{
    int pos = 0;
    if (o->arg == OPTION_OPTIONAL_ARGUMENT || o->arg == OPTION_REQUIRED_ARGUMENT) {
        assert(o->variable != NULL);
        switch (o->type) {
        case ('F'):{
                pos += printf("=%.14f ", *(double *) o->variable);
                break;
            }
        case ('f'):{
                pos += printf("=%.6f ", (double) *(float *) o->variable);
                break;
            }
        case ('d'):{
                pos += printf("=%d ", *(int *) o->variable);
                break;
            }
        case ('H'):
        case ('s'):{
                if (*(char * *) o->variable != NULL && ((char * *) o->variable)[0][0] != 0) {
                    pos += printf("=%s", *(char * *) o->variable);
                } else {
                    pos += printf("=STRING");
                }
                break;
            }
        case ('c'):{
                pos += printf("=%c", *(char *) o->variable);
                break;
            }
        case ('l'):{
                pos += printf("=%lld", *(long long *) o->variable);
                break;
            }
        }
    }
    if (o->arg == OPTION_FLAG && (*(int *) o->variable) != 0) {
        pos += printf(" (%d)", (*(int *) o->variable));
    }
    return pos;
}

static void print_help_section(option_help * args, option_value_type type, char *name)
{
    int first;
    first = 1;
    option_help *o;
    for (o = args; o->shortVar != 0 || o->longVar != 0; o++) {
        if (o->arg == type) {
            if (first) {
                printf("\n%s\n", name);
                first = 0;
            }
            int pos = 0;
            if (o->shortVar != 0 && o->longVar != 0) {
                pos += printf("-%c, --%s", o->shortVar, o->longVar);
            } else if (o->shortVar != 0) {
                pos += printf("-%c", o->shortVar);
            } else if (o->longVar != 0) {
                pos += printf("--%s", o->longVar);
            }
            pos += print_value(o);
            if (o->help != NULL) {
                for (int i = 0; i < (30 - pos); i++) {
                    printf(" ");
                }
                printf("%s", o->help);
            }
            printf("\n");
        }
    }
}

void print_help(option_help * args, int is_plugin)
{
    option_help *o;
    int optionalArgs = 0;
    for (o = args; o->shortVar != 0 || o->longVar != 0; o++) {
        if (o->arg != OPTION_REQUIRED_ARGUMENT) {
            optionalArgs = 1;
        }
        switch (o->arg) {
        case (OPTION_OPTIONAL_ARGUMENT):
        case (OPTION_FLAG):{
                if (o->shortVar != 0) {
                    printf("[-%c] ", o->shortVar);
                } else if (o->longVar != 0) {
                    printf("[--%s] ", o->longVar);
                }
                break;
            }
        case (OPTION_REQUIRED_ARGUMENT):{
                if (o->shortVar != 0) {
                    printf("-%c ", o->shortVar);
                } else if (o->longVar != 0) {
                    printf("--%s ", o->longVar);
                }
                break;
            }
        }
    }
    if (optionalArgs) {
    }
//printf(" [Optional Args]");
    print_help_section(args, OPTION_REQUIRED_ARGUMENT, "Required arguments");
    print_help_section(args, OPTION_FLAG, "Flags");
    print_help_section(args, OPTION_OPTIONAL_ARGUMENT, "Optional arguments");
}

static int print_option_value(option_help * o)
{
    int pos = 0;
    if (o->arg == OPTION_OPTIONAL_ARGUMENT || o->arg == OPTION_REQUIRED_ARGUMENT) {
        assert(o->variable != NULL);
        switch (o->type) {
        case ('F'):{
                pos += printf("=%.14f ", *(double *) o->variable);
                break;
            }
        case ('f'):{
                pos += printf("=%.6f ", (double) *(float *) o->variable);
                break;
            }
        case ('d'):{
                pos += printf("=%d ", *(int *) o->variable);
                break;
            }
        case ('H'):{
                pos += printf("=HIDDEN");
                break;
            }
        case ('s'):{
                if (*(char * *) o->variable != NULL && ((char * *) o->variable)[0][0] != 0) {
                    pos += printf("=%s", *(char * *) o->variable);
                } else {
                    pos += printf("=");
                }
                break;
            }
        case ('c'):{
                pos += printf("=%c", *(char *) o->variable);
                break;
            }
        case ('l'):{
                pos += printf("=%lld", *(long long *) o->variable);
                break;
            }
        }
    } else {
    }
//printf(" ");
    return pos;
}

static void print_current_option_section(option_help * args, option_value_type type)
{
    option_help *o;
    for (o = args; o->shortVar != 0 || o->longVar != 0; o++) {
        if (o->arg == type) {
            int pos = 0;
            if (o->arg == OPTION_FLAG && (*(int *) o->variable) == 0) {
                continue;
            }
            printf("\t");
            if (o->shortVar != 0 && o->longVar != 0) {
                pos += printf("%s", o->longVar);
            } else if (o->shortVar != 0) {
                pos += printf("%c", o->shortVar);
            } else if (o->longVar != 0) {
                pos += printf("%s", o->longVar);
            }
            pos += print_option_value(o);
            printf("\n");
        }
    }
}

void print_current_options(option_help * args)
{
    print_current_option_section(args, OPTION_REQUIRED_ARGUMENT);
    print_current_option_section(args, OPTION_OPTIONAL_ARGUMENT);
    print_current_option_section(args, OPTION_FLAG);
}

int parseOptions(int argc, char * *argv, option_help * args, int *printhelp)
{
    int error = 0;
    int requiredArgsSeen = 0;
    int requiredArgsNeeded = 0;
    int i;
    for (option_help * o = args; o->shortVar != 0 || o->longVar != 0; o++) {
        if (o->arg == OPTION_REQUIRED_ARGUMENT) {
            requiredArgsNeeded++;
        }
    }
    for (i = 1; i < argc; i++) {
        char *txt = argv[i];
        int foundOption = 0;
        char *arg = strstr(txt, "=");
        if (arg != NULL) {
            arg[0] = 0;
            arg++;
        }
        if (strcmp(txt, "--") == 0) {
// we found plugin options
            break;
        }
// try to find matching option help
        for (option_help * o = args; o->shortVar != 0 || o->longVar != 0; o++) {
            if ((strlen(txt) == 2 && txt[0] == '-' && o->shortVar == txt[1]) || (strlen(txt) > 2 && txt[0] == '-' && txt[1] == '-' && o->longVar != NULL && strcmp(txt + 2, o->longVar) == 0)) {
                foundOption = 1;
// now process the option.
                switch (o->arg) {
                case (OPTION_FLAG):{
                        assert(o->type == 'd');
                        (*(int *) o->variable)++;
                        break;
                    }
                case (OPTION_OPTIONAL_ARGUMENT):
                case (OPTION_REQUIRED_ARGUMENT):{
// check if next is an argument
                        if (arg == NULL) {
                            char str_s[1024] = "";
                            if (o->shortVar) {
                                strcat(str_s, "-");
                                strncat(str_s, &o->shortVar, 1);
                                strcat(str_s, " ");
                            }
                            if (o->longVar) {
                                strcat(str_s, "(--");
                                strcat(str_s, o->longVar);
                                strcat(str_s, ") ");
                            }
                            printf("Error option %srequires an argument.\n", str_s);
                            error = 1;
                            break;
                        }
                        switch (o->type) {
                        case ('F'):{
                                *(double *) o->variable = atof(arg);
                                break;
                            }
                        case ('f'):{
                                *(float *) o->variable = atof(arg);
                                break;
                            }
                        case ('d'):{
                                *(int *) o->variable = atoi(arg);
                                break;
                            }
                        case ('H'):
                        case ('s'):{
                                (*(char * *) o->variable) = strdup(arg);
                                break;
                            }
                        case ('c'):{
                                (*(char *) o->variable) = arg[0];
                                if (strlen(arg) > 1) {
                                    printf("Error, ignoring remainder of string for option %c (%s).\n", o->shortVar, o->longVar);
                                }
                                break;
                            }
                        case ('l'):{
                                *(long long *) o->variable = atoll(arg);
                                break;
                            }
                        }
                    }
                }
                if (o->arg == OPTION_REQUIRED_ARGUMENT) {
                    requiredArgsSeen++;
                }
                break;
            }
        }
        if (!foundOption) {
            if (strcmp(txt, "-h") == 0 || strcmp(txt, "--help") == 0) {
                *printhelp = 1;
            } else {
                printf("Error invalid argument: %s\n", txt);
                error = 1;
            }
        }
    }
    if (requiredArgsSeen != requiredArgsNeeded) {
        printf("Error: Missing some required arguments\n\n");
        *printhelp = -1;
    }
    if (error != 0) {
        printf("Invalid options\n");
        *printhelp = -1;
    }
    return i;
}
