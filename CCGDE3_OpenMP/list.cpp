/*
  list.c

  Author:
       Luis Miguel Antonio <lmiguel@computacion.cs.cinvestav.mx>       

  Copyright (c) 2013 Luis Miguel Antonio		*/
  
  
#include "CCGDE3.h"

void CCGDE3::insert(node *n, int x){
    node *temp;
    if (n==NULL){
        printf("\n Error!! NULL pointer\n");
        exit(1);
    }
    temp = (node *)calloc(1,sizeof(node));
    temp->index = x;
    temp->child = n->child;
    temp->parent = n;
    if (n->child != NULL){
        n->child->parent = temp;
    }
    n->child = temp;
    return;
}

list* CCGDE3::createList(int index){
	list *l;
	l = (list *)calloc(1,sizeof(list));
	l->index = index;
	l->parent = NULL;
	l->child = NULL;
	return l;
}


node* CCGDE3::deleteNode(node *n){
    node *temp;
    if(n==NULL){
        printf("\n Error!! puntero NULL\n");
        exit(1);
    }
    temp = n->parent;
    temp->child = n->child;
    if (temp->child!=NULL){
        temp->child->parent = temp;
    }
    free (n);
    return (temp);
}


node* CCGDE3::deleteInd(list *l,int index){
    node *temp;
    if (l==NULL){
        printf("\n Error!! NULL pointer\n");
        exit(1);
    }
    temp = l->child;
    while(temp != NULL){
		if(temp->index==index){
			break;
		}
		temp = temp->child;
	}
	if(temp == NULL){
		printf("\n Error!! no index in list\n");
        exit(1);
	}
	temp = deleteNode(temp);
    return (temp);
}

void CCGDE3::deleteList(list *l){
	node *temp1;
	while (l != NULL){
        temp1 = l;
        l = l->child;
        free(temp1);
    }
}
